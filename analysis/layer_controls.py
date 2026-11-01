#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""B1 + B2：三层对照与"仅扩增支撑"拆分（V3 补充分析）

B1 三层对照（同一套 30 基因判定）
  L1 原发 vs 转移（跨界）           —— 主分析
  L2 转移 vs 转移，同器官（噪声地板）
  L3 转移 vs 转移，同患者不同器官（新增，检验"器官身份"本身的作用）
  说明：L3 在每个患者内只取一对（按器官字母序 + 覆盖度），保证患者间独立。

B2 拆分
  (a) 仅由扩增支撑的"漏检"（met-only / primary-only）
  (b) SNV/indel-only 定义下的重复计算（去扩增），含标准化

输出：tables/layer_controls.csv, tables/amplification_decomposition.csv
"""
import csv, collections, os
import numpy as np

BASE = r"F:\数据库汇总\msk_impact_50k_2026"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "05_results", "tables")
SEED = 20261101

E2S = {"1956":"EGFR","3845":"KRAS","673":"BRAF","2064":"ERBB2","4233":"MET","5290":"PIK3CA","672":"BRCA1",
       "675":"BRCA2","3417":"IDH1","3418":"IDH2","2260":"FGFR1","2263":"FGFR2","2261":"FGFR3","5727":"PTCH1",
       "6608":"SMO","207":"AKT1","2475":"MTOR","7248":"TSC1","7249":"TSC2","3815":"KIT","5156":"PDGFRA",
       "25":"ABL1","5979":"RET","238":"ALK","6098":"ROS1","4914":"NTRK1","2099":"ESR1","367":"AR",
       "1019":"CDK4","4193":"MDM2"}
GENES = set(E2S.values())
NONCODING = {"Silent","Intron","3'UTR","5'UTR","IGR","RNA","Flank","5'Flank","3'Flank"}

smp = {r["sampleId"]: r for r in csv.DictReader(open(os.path.join(BASE, "data_clinical_sample.csv"), encoding="utf-8-sig"))}
pid = lambda s: "-".join(s.split("-")[:2])

snv = collections.defaultdict(set)
for r in csv.DictReader(open(os.path.join(BASE, "data_mutations.csv"), encoding="utf-8-sig")):
    g = E2S.get(r["entrezGeneId"])
    if g and r["mutationType"] not in NONCODING:
        snv[r["sampleId"]].add(g)

amp = collections.defaultdict(set)
with open(os.path.join(BASE, "data_cna.txt"), encoding="utf-8-sig") as f:
    head = f.readline().rstrip("\n").split("\t")[1:]
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[0] not in GENES:
            continue
        for sid, v in zip(head, parts[1:]):
            try:
                if float(v) >= 2:
                    amp[sid].add(parts[0])
            except ValueError:
                pass

any_alt = lambda s: bool(snv.get(s) or amp.get(s))
snv_only = lambda s: bool(snv.get(s))

def cov(x):
    try: return int(x["SAMPLE_COVERAGE"] or 0)
    except: return 0

# ---------- 队列
bypat = collections.defaultdict(list)
for s, r in smp.items():
    bypat[pid(s)].append(r)

pairs = []
for p, rs in bypat.items():
    prim = [x for x in rs if x["SAMPLE_TYPE"] == "Primary"]
    met = [x for x in rs if x["SAMPLE_TYPE"] == "Metastasis"]
    if prim and met:
        prim.sort(key=cov, reverse=True); met.sort(key=cov, reverse=True)
        pairs.append((p, prim[0], met[0]))

same_organ = []
for p, rs in bypat.items():
    d = collections.defaultdict(list)
    for x in rs:
        if x["SAMPLE_TYPE"] == "Metastasis" and x["METASTATIC_SITE"].strip():
            d[x["METASTATIC_SITE"].strip().lower()].append(x)
    for organ, v in d.items():
        if len(v) >= 2:
            v.sort(key=cov, reverse=True)
            same_organ.append((p, organ, v[0], v[1]))

# L3：同患者不同器官，每患者 1 对
cross_organ = []
for p, rs in bypat.items():
    d = collections.defaultdict(list)
    for x in rs:
        if x["SAMPLE_TYPE"] == "Metastasis" and x["METASTATIC_SITE"].strip():
            d[x["METASTATIC_SITE"].strip().lower()].append(x)
    if len(d) >= 2:
        orgs = sorted(d)
        a = sorted(d[orgs[0]], key=cov, reverse=True)[0]
        b = sorted(d[orgs[1]], key=cov, reverse=True)[0]
        cross_organ.append((p, orgs[0], orgs[1], a, b))

def rate(flags):
    k = sum(flags); n = len(flags)
    return k, n, (100 * k / n if n else float("nan"))

def std_rate(flags_by_cancer, weights_from):
    """按 weights_from（癌种->权重）标准化"""
    num = den = 0.0
    for c, w in weights_from.items():
        fl = flags_by_cancer.get(c)
        if fl is None or len(fl) < 5:
            continue
        num += w * np.mean(fl); den += w
    return 100 * num / den if den else float("nan")

# ---------- L1/L2/L3（any alteration）
l1 = [int(any_alt(a["sampleId"]) != any_alt(b["sampleId"])) for _, a, b in pairs]
l2 = [int(any_alt(a["sampleId"]) != any_alt(b["sampleId"])) for _, _, a, b in same_organ]
l3 = [int(any_alt(a["sampleId"]) != any_alt(b["sampleId"])) for _, _, _, a, b in cross_organ]

# 标准化权重 = L1 的癌种分布
w = collections.Counter(a["CANCER_TYPE"] for _, a, b in pairs)
w = {c: v / sum(w.values()) for c, v in w.items()}
by_c_L1 = collections.defaultdict(list)
for _, a, b in pairs: by_c_L1[a["CANCER_TYPE"]].append(int(any_alt(a["sampleId"]) != any_alt(b["sampleId"])))
by_c_L2 = collections.defaultdict(list)
for _, _, a, b in same_organ: by_c_L2[smp[a["sampleId"]]["CANCER_TYPE"]].append(int(any_alt(a["sampleId"]) != any_alt(b["sampleId"])))
by_c_L3 = collections.defaultdict(list)
for _, _, _, a, b in cross_organ: by_c_L3[smp[a["sampleId"]]["CANCER_TYPE"]].append(int(any_alt(a["sampleId"]) != any_alt(b["sampleId"])))

rows = []
for name, flags, byc in [("L1_primary_vs_metastasis", l1, by_c_L1),
                         ("L2_same_organ_replicates", l2, by_c_L2),
                         ("L3_same_patient_different_organs", l3, by_c_L3)]:
    k, n, pct = rate(flags)
    rows.append({"layer": name, "n": n, "discordant": k, "pct": round(pct, 2),
                 "pct_standardised": round(std_rate(byc, w), 2)})

# SNV-only（去扩增）重复
l1s = [int(snv_only(a["sampleId"]) != snv_only(b["sampleId"])) for _, a, b in pairs]
l2s = [int(snv_only(a["sampleId"]) != snv_only(b["sampleId"])) for _, _, a, b in same_organ]
l3s = [int(snv_only(a["sampleId"]) != snv_only(b["sampleId"])) for _, _, _, a, b in cross_organ]
by1 = collections.defaultdict(list); by2 = collections.defaultdict(list); by3 = collections.defaultdict(list)
for _, a, b in pairs: by1[a["CANCER_TYPE"]].append(int(snv_only(a["sampleId"]) != snv_only(b["sampleId"])))
for _, _, a, b in same_organ: by2[smp[a["sampleId"]]["CANCER_TYPE"]].append(int(snv_only(a["sampleId"]) != snv_only(b["sampleId"])))
for _, _, _, a, b in cross_organ: by3[smp[a["sampleId"]]["CANCER_TYPE"]].append(int(snv_only(a["sampleId"]) != snv_only(b["sampleId"])))
for name, flags, byc in [("L1_SNVonly", l1s, by1), ("L2_SNVonly", l2s, by2), ("L3_SNVonly", l3s, by3)]:
    k, n, pct = rate(flags)
    rows.append({"layer": name, "n": n, "discordant": k, "pct": round(pct, 2),
                 "pct_standardised": round(std_rate(byc, w), 2)})

with open(os.path.join(TAB, "layer_controls.csv"), "w", newline="", encoding="utf-8-sig") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)

# ---------- B2 拆分
met_only = [(p, a, b) for p, a, b in pairs if any_alt(b["sampleId"]) and not any_alt(a["sampleId"])]
prim_only = [(p, a, b) for p, a, b in pairs if any_alt(a["sampleId"]) and not any_alt(b["sampleId"])]

def classify(sel, which):
    amp_only = snv_and_amp = snv_only_n = 0
    for p, a, b in sel:
        s = b["sampleId"] if which == "met" else a["sampleId"]
        has_s, has_a = bool(snv.get(s)), bool(amp.get(s))
        if has_a and not has_s: amp_only += 1
        elif has_a and has_s: snv_and_amp += 1
        else: snv_only_n += 1
    return amp_only, snv_and_amp, snv_only_n

m = classify(met_only, "met"); pr = classify(prim_only, "prim")
decomp = [
    {"direction": "met_only (lost if primary-tested)", "total": len(met_only),
     "amplification_only": m[0], "snv_and_amplification": m[1], "snv_only": m[2],
     "amp_only_pct": round(100 * m[0] / len(met_only), 1)},
    {"direction": "primary_only (lost if metastasis-tested)", "total": len(prim_only),
     "amplification_only": pr[0], "snv_and_amplification": pr[1], "snv_only": pr[2],
     "amp_only_pct": round(100 * pr[0] / len(prim_only), 1)},
]
with open(os.path.join(TAB, "amplification_decomposition.csv"), "w", newline="", encoding="utf-8-sig") as f:
    wr = csv.DictWriter(f, fieldnames=list(decomp[0].keys())); wr.writeheader(); wr.writerows(decomp)

print("=== layer controls ===")
for r in rows: print(f"  {r['layer']:38s} n={r['n']:5d} disc={r['discordant']:4d} {r['pct']:5.1f}% | standardised {r['pct_standardised']:5.1f}%")
print("=== amplification decomposition ===")
for d in decomp: print(f"  {d['direction']:42s} total={d['total']:4d} amp-only={d['amplification_only']:3d} ({d['amp_only_pct']}%) snv-only={d['snv_only']:3d} both={d['snv_and_amplification']:3d}")
