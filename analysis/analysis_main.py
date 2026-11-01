#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
OIR 项目 主分析脚本（analysis_plan_v2.0 对应实现）
输入：F:\\数据库汇总\\msk_impact_50k_2026
输出：05_results/tables/*.csv, 05_results/results_summary.json, 05_results/statistical_report.md

设计要点（与 v2.0 一一对应）：
  Step 1  队列构建：配对对（原发 vs 转移）、同器官重复对（噪声地板）
  Step 2  可靶向状态：30 基因 SNV/indel（排除同义等）+ 扩增（CNA>=2）
  Step 3  P1 总体翻转率 + 噪声地板 + 器官级超额翻转（bootstrap CI）
  Step 4  P2 决策影响（仅测原发/仅测转移会漏掉什么）
  Step 5  TMB-H / MSI-H 一致性（含噪声地板）
  Step 6  P5 OS（Cox，按癌种分层；证伪性检验）
随机种子固定：20261101
"""
from __future__ import annotations
import csv, json, os, random, statistics, sys
from collections import defaultdict, Counter
from datetime import datetime

import numpy as np
import pandas as pd

SEED = 20261101
random.seed(SEED); np.random.seed(SEED)

BASE = r"F:\数据库汇总\msk_impact_50k_2026"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "05_results")
TAB = os.path.join(OUT, "tables")
os.makedirs(TAB, exist_ok=True)

TARGET_GENES = {
    "EGFR","KRAS","BRAF","ERBB2","MET","PIK3CA","BRCA1","BRCA2","IDH1","IDH2",
    "FGFR1","FGFR2","FGFR3","PTCH1","SMO","AKT1","MTOR","TSC1","TSC2","KIT",
    "PDGFRA","ABL1","RET","ALK","ROS1","NTRK1","ESR1","AR","CDK4","MDM2",
}
# entrez -> symbol（仅本分析所需）
E2S = {
    "1956":"EGFR","3845":"KRAS","673":"BRAF","2064":"ERBB2","4233":"MET","5290":"PIK3CA",
    "672":"BRCA1","675":"BRCA2","3417":"IDH1","3418":"IDH2","2260":"FGFR1","2263":"FGFR2",
    "2261":"FGFR3","5727":"PTCH1","6608":"SMO","207":"AKT1","2475":"MTOR","7248":"TSC1",
    "7249":"TSC2","3815":"KIT","5156":"PDGFRA","25":"ABL1","5979":"RET","238":"ALK",
    "6098":"ROS1","4914":"NTRK1","2099":"ESR1","367":"AR","1019":"CDK4","4193":"MDM2",
}
NONCODING = {"Silent","Intron","3'UTR","5'UTR","IGR","RNA","Flank","5'Flank","3'Flank"}

def log(msg): print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)

# ---------------------------------------------------------------- Step 1 队列
def load_clinical():
    samples = {}
    with open(os.path.join(BASE, "data_clinical_sample.csv"), encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            samples[r["sampleId"]] = r
    patients = {}
    with open(os.path.join(BASE, "data_clinical_patient.txt"), encoding="utf-8-sig") as f:
        hdr = None
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "PATIENT_ID": hdr = parts; continue
            patients[parts[0]] = dict(zip(hdr, parts))
    return samples, patients

def pid_of(sample_id): return "-".join(sample_id.split("-")[:2])

def cov(s):
    try: return int(s.get("SAMPLE_COVERAGE") or 0)
    except: return 0

def build_pairs(samples):
    bypat = defaultdict(list)
    for sid, r in samples.items():
        bypat[pid_of(sid)].append(r)
    pairs = []
    for p, rs in bypat.items():
        prim = [x for x in rs if x["SAMPLE_TYPE"] == "Primary"]
        met  = [x for x in rs if x["SAMPLE_TYPE"] == "Metastasis"]
        if prim and met:
            prim.sort(key=cov, reverse=True); met.sort(key=cov, reverse=True)
            pairs.append({"patient": p, "primary": prim[0], "met": met[0],
                          "organ": (met[0]["METASTATIC_SITE"] or "Unspecified").strip(),
                          "cancer": prim[0]["CANCER_TYPE"]})
    return pairs

def build_noise_pairs(samples):
    """同器官重复对：同一患者、同一 METASTATIC_SITE、>=2 个转移样本，取覆盖度前二"""
    by = defaultdict(list)
    for sid, r in samples.items():
        if r["SAMPLE_TYPE"] == "Metastasis" and r["METASTATIC_SITE"].strip():
            by[(pid_of(sid), r["METASTATIC_SITE"].strip().lower())].append(r)
    out = []
    for (p, organ), v in by.items():
        if len(v) >= 2:
            v.sort(key=cov, reverse=True)
            out.append({"patient": p, "organ": organ.title(), "a": v[0], "b": v[1]})
    return out

# ------------------------------------------------------- Step 2 可靶向状态
def load_alterations():
    """sample -> set(gene symbols) 的可靶向改变（SNV/indel + 扩增）"""
    alt = defaultdict(set)
    n_rows = 0
    with open(os.path.join(BASE, "data_mutations.csv"), encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            n_rows += 1
            g = E2S.get(r["entrezGeneId"])
            if g and r["mutationType"] not in NONCODING:
                alt[r["sampleId"]].add(g)
    log(f"mutations rows={n_rows}; samples with SNV/indel alteration={len(alt)}")

    # CNA：只保留 30 基因所在行，值 >=2 记为扩增
    cna_path = os.path.join(BASE, "data_cna.txt")
    n_amp = 0
    with open(cna_path, encoding="utf-8-sig") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        sample_ids = hdr[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]
            if gene not in TARGET_GENES: continue
            for sid, val in zip(sample_ids, parts[1:]):
                if val in ("", "NA", "nan"): continue
                try:
                    if float(val) >= 2:
                        alt[sid].add(gene); n_amp += 1
                except ValueError:
                    continue
    log(f"CNA amplification calls (>=2)={n_amp}")
    return alt

# ---------------------------------------------------------------- 工具函数
def wilson(k, n, z=1.96):
    if n == 0: return (float("nan"),) * 3
    p = k / n
    d = 1 + z*z/n
    c = p + z*z/(2*n)
    m = z*((p*(1-p)/n + z*z/(4*n*n)) ** 0.5)
    return (p*100, (c-m)/d*100, (c+m)/d*100)

def boot_excess(cross_flags, noise_flags, n_boot=2000):
    """超额 = P(cross) - P(noise)，患者级 bootstrap"""
    if not cross_flags or not noise_flags: return (float("nan"),)*3
    c = np.array(cross_flags, dtype=float); n = np.array(noise_flags, dtype=float)
    obs = 100*(c.mean() - n.mean())
    diffs = []
    rng = np.random.default_rng(SEED)
    for _ in range(n_boot):
        ci = rng.integers(0, len(c), len(c)); ni = rng.integers(0, len(n), len(n))
        diffs.append(100*(c[ci].mean() - n[ni].mean()))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return obs, lo, hi

def tmb_of(s):
    try: return float(s["TMB_SCORE"])
    except: return None

def is_tmb_high(s, cut=10.0):
    v = tmb_of(s)
    return None if v is None else (v >= cut)

def main():
    log("Step 1 载入临床数据")
    samples, patients = load_clinical()
    log(f"samples={len(samples)} patients={len(patients)}")
    pairs = build_pairs(samples)
    noise = build_noise_pairs(samples)
    log(f"pairs={len(pairs)}  same-organ replicate pairs={len(noise)}")

    log("Step 2 载入可靶向改变（SNV/indel + 扩增）")
    alt = load_alterations()

    # ---- P1 总体
    cross_flags = []
    for x in pairs:
        a = x["primary"]["sampleId"] in alt
        b = x["met"]["sampleId"] in alt
        cross_flags.append(a != b)
    noise_flags = []
    for x in noise:
        a = x["a"]["sampleId"] in alt
        b = x["b"]["sampleId"] in alt
        noise_flags.append(a != b)

    p_cross = wilson(sum(cross_flags), len(cross_flags))
    p_noise = wilson(sum(noise_flags), len(noise_flags))
    excess = boot_excess(cross_flags, noise_flags)
    log(f"P1 cross={sum(cross_flags)}/{len(cross_flags)} {p_cross[0]:.1f}%  "
        f"noise={sum(noise_flags)}/{len(noise_flags)} {p_noise[0]:.1f}%  "
        f"excess={excess[0]:.1f}pp ({excess[1]:.1f}–{excess[2]:.1f})")

    rows = [{
        "metric": "gene_level_discordance",
        "cross_n": len(cross_flags), "cross_k": sum(cross_flags),
        "cross_pct": round(p_cross[0], 2), "cross_lo": round(p_cross[1], 2), "cross_hi": round(p_cross[2], 2),
        "noise_n": len(noise_flags), "noise_k": sum(noise_flags),
        "noise_pct": round(p_noise[0], 2), "noise_lo": round(p_noise[1], 2), "noise_hi": round(p_noise[2], 2),
        "excess_pp": round(excess[0], 2), "excess_lo": round(excess[1], 2), "excess_hi": round(excess[2], 2),
    }]

    # ---- P1 器官级
    cross_by_org = defaultdict(list); noise_by_org = defaultdict(list)
    for x in pairs:
        cross_by_org[x["organ"].title()].append(
            (x["primary"]["sampleId"] in alt) != (x["met"]["sampleId"] in alt))
    for x in noise:
        noise_by_org[x["organ"].title()].append(
            (x["a"]["sampleId"] in alt) != (x["b"]["sampleId"] in alt))

    organ_rows = []
    for org, flags in sorted(cross_by_org.items(), key=lambda kv: -len(kv[1])):
        nf = noise_by_org.get(org, [])
        if len(flags) < 30: continue
        c = wilson(sum(flags), len(flags))
        nz = wilson(sum(nf), len(nf)) if nf else (float("nan"),)*3
        ex = boot_excess(flags, nf) if len(nf) >= 15 else (float("nan"),)*3
        organ_rows.append({
            "organ": org, "cross_n": len(flags), "cross_k": sum(flags), "cross_pct": round(c[0], 2),
            "cross_lo": round(c[1], 2), "cross_hi": round(c[2], 2),
            "noise_n": len(nf), "noise_k": sum(nf),
            "noise_pct": round(nz[0], 2) if nf else None,
            "excess_pp": round(ex[0], 2) if not np.isnan(ex[0]) else None,
            "excess_lo": round(ex[1], 2) if not np.isnan(ex[1]) else None,
            "excess_hi": round(ex[2], 2) if not np.isnan(ex[2]) else None,
        })
    pd.DataFrame(organ_rows).to_csv(os.path.join(TAB, "p1_organ_excess_discordance.csv"), index=False)
    log(f"P1 器官级：{len(organ_rows)} 个器官（cross n>=30）")

    # ---- 敏感性分析（同 panel；双端覆盖度>=200×）
    sens_rows = []
    def subset_stats(sel_pairs, sel_noise, label):
        cf = [ (x["primary"]["sampleId"] in alt) != (x["met"]["sampleId"] in alt) for x in sel_pairs ]
        nf = [ (x["a"]["sampleId"] in alt) != (x["b"]["sampleId"] in alt) for x in sel_noise ]
        if not cf or not nf:
            sens_rows.append({"subset": label, "cross_n": len(cf), "cross_pct": None,
                              "noise_n": len(nf), "noise_pct": None, "excess_pp": None}); return
        c = wilson(sum(cf), len(cf)); nz = wilson(sum(nf), len(nf)); ex = boot_excess(cf, nf)
        sens_rows.append({"subset": label,
                          "cross_n": len(cf), "cross_k": sum(cf), "cross_pct": round(c[0], 2),
                          "noise_n": len(nf), "noise_k": sum(nf), "noise_pct": round(nz[0], 2),
                          "excess_pp": round(ex[0], 2), "excess_lo": round(ex[1], 2), "excess_hi": round(ex[2], 2)})
        log(f"  敏感性 {label}: cross {c[0]:.1f}% noise {nz[0]:.1f}% excess {ex[0]:.1f}pp ({ex[1]:.1f}–{ex[2]:.1f})")

    subset_stats([x for x in pairs if x["primary"]["GENE_PANEL"] == x["met"]["GENE_PANEL"]],
                 [x for x in noise if x["a"]["GENE_PANEL"] == x["b"]["GENE_PANEL"]], "same_panel")
    subset_stats([x for x in pairs if cov(x["primary"]) >= 200 and cov(x["met"]) >= 200],
                 [x for x in noise if cov(x["a"]) >= 200 and cov(x["b"]) >= 200], "coverage_ge200")
    pd.DataFrame(sens_rows).to_csv(os.path.join(TAB, "p1_sensitivity.csv"), index=False)

    # ---- 器官可替代性规则（超额 CI 上限 <5pp 且 cross n>=50 → 可替代）
    rule_rows = []
    for r in organ_rows:
        if r["excess_pp"] is None or r["cross_n"] < 50:
            verdict = "insufficient"
        elif r["excess_hi"] is not None and r["excess_hi"] < 5:
            verdict = "substitutable"
        elif r["excess_lo"] is not None and r["excess_lo"] > 5:
            verdict = "not_substitutable"
        else:
            verdict = "uncertain"
        rule_rows.append({**r, "rule_verdict": verdict})
    pd.DataFrame(rule_rows).to_csv(os.path.join(TAB, "p4_organ_rule.csv"), index=False)
    for r in rule_rows:
        log(f"  规则 {r['organ']}: {r['rule_verdict']}")

    # ---- P2 决策影响
    met_only = [x for x in pairs if (x["met"]["sampleId"] in alt) and not (x["primary"]["sampleId"] in alt)]
    prim_only = [x for x in pairs if (x["primary"]["sampleId"] in alt) and not (x["met"]["sampleId"] in alt)]
    gene_gain, gene_loss = Counter(), Counter()
    for x in met_only:
        for g in alt[x["met"]["sampleId"]]: gene_gain[g] += 1
    for x in prim_only:
        for g in alt[x["primary"]["sampleId"]]: gene_loss[g] += 1

    p2 = [{
        "metric": "missed_if_primary_only", "k": len(met_only), "n": len(pairs),
        "pct": round(wilson(len(met_only), len(pairs))[0], 2),
        "lo": round(wilson(len(met_only), len(pairs))[1], 2),
        "hi": round(wilson(len(met_only), len(pairs))[2], 2),
    }, {
        "metric": "missed_if_metastasis_only", "k": len(prim_only), "n": len(pairs),
        "pct": round(wilson(len(prim_only), len(pairs))[0], 2),
        "lo": round(wilson(len(prim_only), len(pairs))[1], 2),
        "hi": round(wilson(len(prim_only), len(pairs))[2], 2),
    }]
    pd.DataFrame(p2).to_csv(os.path.join(TAB, "p2_decision_impact.csv"), index=False)
    pd.DataFrame([{"gene": g, "patients_missed_if_primary_only": n} for g, n in gene_gain.most_common()]
                 ).to_csv(os.path.join(TAB, "p2_gene_gain_in_metastasis.csv"), index=False)
    pd.DataFrame([{"gene": g, "patients_missed_if_metastasis_only": n} for g, n in gene_loss.most_common()]
                 ).to_csv(os.path.join(TAB, "p2_gene_lost_in_metastasis.csv"), index=False)
    log(f"P2 met-only={len(met_only)} ({p2[0]['pct']}%)  primary-only={len(prim_only)} ({p2[1]['pct']}%)")

    # ---- Step 5 TMB-H / MSI-H
    tmb_rows = []
    t_pairs = [(is_tmb_high(x["primary"]), is_tmb_high(x["met"])) for x in pairs]
    t_pairs = [(a, b) for a, b in t_pairs if a is not None and b is not None]
    k = sum(1 for a, b in t_pairs if a != b)
    c = wilson(k, len(t_pairs))
    t_noise = [(is_tmb_high(x["a"]), is_tmb_high(x["b"])) for x in noise]
    t_noise = [(a, b) for a, b in t_noise if a is not None and b is not None]
    kn = sum(1 for a, b in t_noise if a != b)
    cn = wilson(kn, len(t_noise))
    ex = boot_excess([a != b for a, b in t_pairs], [a != b for a, b in t_noise])
    tmb_rows.append({"metric": "TMB_H_10", "cross_n": len(t_pairs), "cross_k": k,
                     "cross_pct": round(c[0], 2), "cross_lo": round(c[1], 2), "cross_hi": round(c[2], 2),
                     "noise_n": len(t_noise), "noise_k": kn, "noise_pct": round(cn[0], 2),
                     "noise_lo": round(cn[1], 2), "noise_hi": round(cn[2], 2),
                     "excess_pp": round(ex[0], 2), "excess_lo": round(ex[1], 2), "excess_hi": round(ex[2], 2)})
    pd.DataFrame(tmb_rows).to_csv(os.path.join(TAB, "tmb_discordance.csv"), index=False)

    # MSI
    msi_pairs = [(x["primary"]["MSI_TYPE"], x["met"]["MSI_TYPE"]) for x in pairs]
    msi_pairs = [(a, b) for a, b in msi_pairs if a in ("Stable", "Instable") and b in ("Stable", "Instable")]
    km = sum(1 for a, b in msi_pairs if (a == "Instable") != (b == "Instable"))
    pd.DataFrame([{"metric": "MSI_H", "n": len(msi_pairs), "k": km,
                   "pct": round(wilson(km, len(msi_pairs))[0], 2)}]
                 ).to_csv(os.path.join(TAB, "msi_discordance.csv"), index=False)
    log(f"TMB-H cross={k}/{len(t_pairs)} noise={kn}/{len(t_noise)}; MSI-H={km}/{len(msi_pairs)}")

    # ---- Step 6 OS（Cox）
    os_rows = []
    try:
        from lifelines import CoxPHFitter
        from lifelines.statistics import proportional_hazard_test
        recs = []
        for x in pairs:
            p = x["patient"]; pr = patients.get(p, {})
            try:
                t = float(pr["OS_MONTHS"])
            except (KeyError, ValueError):
                continue
            ev = 1 if pr.get("OS_STATUS", "").strip().upper() == "DECEASED" else 0
            tp, tm = tmb_of(x["primary"]), tmb_of(x["met"])
            tmb_mean = np.mean([v for v in (tp, tm) if v is not None]) if (tp or tm) else np.nan
            try: fga = float(x["met"]["FRACTION_GENOME_ALTERED"])
            except (KeyError, ValueError): fga = np.nan
            try: age = float(pr.get("AGE_AT_DX", "nan"))
            except ValueError: age = np.nan
            recs.append({
                "duration": t, "event": ev,
                "discordant": int((x["primary"]["sampleId"] in alt) != (x["met"]["sampleId"] in alt)),
                "cancer": x["cancer"], "sex": 1 if pr.get("SEX", "").upper().startswith("M") else 0,
                "age": age, "tmb": tmb_mean, "fga": fga,
                "same_panel": int(x["primary"]["GENE_PANEL"] == x["met"]["GENE_PANEL"]),
            })
        df = pd.DataFrame(recs).dropna(subset=["duration", "event", "discordant"])
        df["event"] = df["event"].astype(int)
        log(f"OS 分析集 n={len(df)} events={int(df['event'].sum())}")

        # 描述性中位 OS
        for lab, sub in [("discordant", df[df.discordant == 1]), ("concordant", df[df.discordant == 0])]:
            os_rows.append({"analysis": "descriptive", "group": lab, "n": len(sub),
                            "median_os_months": round(float(np.median(sub["duration"])), 1),
                            "deaths": int(sub["event"].sum()),
                            "death_pct": round(sub["event"].mean() * 100, 1)})

        cph_df = df.copy()
        for col in ["age", "tmb", "fga"]:
            cph_df[col] = cph_df[col].fillna(cph_df[col].median())
        # 癌种按分层处理（避免 66 个哑变量导致完全分离）；样本量<20 的癌种并入 Other
        counts = cph_df["cancer"].value_counts()
        keep = set(counts[counts >= 20].index)
        cph_df["cancer_grp"] = cph_df["cancer"].where(cph_df["cancer"].isin(keep), "Other")
        log(f"  Cox 分层：{cph_df['cancer_grp'].nunique()} 层（>=20 例的癌种 {len(keep)} 个 + Other）")
        cph = CoxPHFitter()
        cph.fit(cph_df, duration_col="duration", event_col="event", strata=["cancer_grp"],
                formula="discordant + age + sex + tmb + fga + same_panel", show_progress=False)
        s = cph.summary.loc["discordant"]
        os_rows.append({"analysis": "cox_adjusted", "group": "discordant_vs_concordant",
                        "n": len(cph_df), "HR": round(float(s["exp(coef)"]), 3),
                        "HR_lo": round(float(s["exp(coef) lower 95%"]), 3),
                        "HR_hi": round(float(s["exp(coef) upper 95%"]), 3),
                        "p": float(s["p"])})
        try:
            ph = proportional_hazard_test(cph, cph_df, time_transform="rank")
            os_rows.append({"analysis": "schoenfeld", "group": "discordant",
                            "p": float(ph.summary.loc["discordant", "p"])})
        except Exception as e:
            log(f"Schoenfeld 检验失败：{e}")
    except ImportError:
        log("lifelines 不可用，跳过 OS")
    pd.DataFrame(os_rows).to_csv(os.path.join(TAB, "os_cox.csv"), index=False)

    # ---- 汇总
    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "seed": SEED,
        "n_pairs": len(pairs), "n_noise_pairs": len(noise),
        "P1_gene_level": rows[0],
        "P2": p2,
        "TMB": tmb_rows[0],
        "MSI": {"n": len(msi_pairs), "k": km},
        "OS": os_rows,
    }
    with open(os.path.join(OUT, "results_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    log("写出 results_summary.json")

    # 简要报告
    with open(os.path.join(OUT, "statistical_report.md"), "w", encoding="utf-8") as f:
        f.write("# 统计结果报告（自动生成，勿手改）\n\n")
        f.write(f"- 生成时间：{summary['generated']}\n- 随机种子：{SEED}\n")
        f.write(f"- 配对：{len(pairs)}；同器官重复对：{len(noise)}\n\n")
        f.write("## P1 主终点\n\n")
        f.write(f"- 跨界翻转：{rows[0]['cross_k']}/{rows[0]['cross_n']} = {rows[0]['cross_pct']}% "
                f"(95%CI {rows[0]['cross_lo']}–{rows[0]['cross_hi']})\n")
        f.write(f"- 噪声地板：{rows[0]['noise_k']}/{rows[0]['noise_n']} = {rows[0]['noise_pct']}% "
                f"(95%CI {rows[0]['noise_lo']}–{rows[0]['noise_hi']})\n")
        f.write(f"- **超额：{rows[0]['excess_pp']} pp (bootstrap 95%CI {rows[0]['excess_lo']}–{rows[0]['excess_hi']})**\n\n")
        f.write("## P2 决策影响\n\n")
        for r in p2:
            f.write(f"- {r['metric']}: {r['k']}/{r['n']} = {r['pct']}% (95%CI {r['lo']}–{r['hi']})\n")
        f.write("\n## TMB-H\n\n")
        f.write(f"- 跨界 {k}/{len(t_pairs)} = {c[0]:.1f}%；噪声地板 {kn}/{len(t_noise)} = {cn[0]:.1f}%；"
                f"超额 {ex[0]:.1f} pp (95%CI {ex[1]:.1f}–{ex[2]:.1f})\n")
        f.write("\n## OS\n\n")
        for r in os_rows:
            f.write(f"- {r}\n")
    log("完成")

if __name__ == "__main__":
    main()
