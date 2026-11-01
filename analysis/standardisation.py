#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A1：癌种标准化超额（V2 主终点）

背景：V1 的主终点直接比较"原发-转移"与"同器官重复"两组，但两组来自不同患者子集，
      癌种构成不同（如前列腺癌 38.9% vs 17.1%）。V2 改为按索引队列癌种分布直接标准化。

方法：对同时具备两类数据的癌种（噪声对数 >=5），以索引队列该癌种占比为权重，
      分别计算加权 cross 率与加权 noise 率，差值即标准化超额；患者级 bootstrap（2000 次）取 CI。
输出：05_results/tables/p1_standardised_excess.csv
"""
import csv, collections, os
import numpy as np

BASE = r"F:\数据库汇总\msk_impact_50k_2026"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "05_results", "tables")
SEED = 20261101

TARGET = {"1956","3845","673","2064","4233","5290","672","675","3417","3418","2260","2263","2261",
          "5727","6608","207","2475","7248","7249","3815","5156","25","5979","238","6098","4914",
          "2099","367","1019","4193"}
NONCODING = {"Silent","Intron","3'UTR","5'UTR","IGR","RNA","Flank","5'Flank","3'Flank"}
CNA_GENES = {"EGFR","ERBB2","MET","MDM2","CDK4","FGFR1","FGFR2","FGFR3","KIT","PIK3CA","AR","ESR1"}


def load():
    smp = {r["sampleId"]: r for r in csv.DictReader(open(os.path.join(BASE, "data_clinical_sample.csv"), encoding="utf-8-sig"))}
    alt = collections.defaultdict(set)
    for r in csv.DictReader(open(os.path.join(BASE, "data_mutations.csv"), encoding="utf-8-sig")):
        if r["entrezGeneId"] in TARGET and r["mutationType"] not in NONCODING:
            alt[r["sampleId"]].add(r["entrezGeneId"])
    with open(os.path.join(BASE, "data_cna.txt"), encoding="utf-8-sig") as f:
        head = f.readline().rstrip("\n").split("\t")[1:]
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] not in CNA_GENES:
                continue
            for sid, v in zip(head, p[1:]):
                try:
                    if float(v) >= 2:
                        alt[sid].add(p[0])
                except ValueError:
                    pass
    return smp, alt


def main():
    smp, alt = load()
    pid = lambda s: "-".join(s.split("-")[:2])
    byp = collections.defaultdict(list)
    for s, r in smp.items():
        byp[pid(s)].append(r)

    pairs = []
    for p, rs in byp.items():
        prim = [x for x in rs if x["SAMPLE_TYPE"] == "Primary"]
        met = [x for x in rs if x["SAMPLE_TYPE"] == "Metastasis"]
        if prim and met:
            prim.sort(key=lambda x: int(x["SAMPLE_COVERAGE"] or 0), reverse=True)
            met.sort(key=lambda x: int(x["SAMPLE_COVERAGE"] or 0), reverse=True)
            pairs.append((p, prim[0], met[0]))

    # 同器官重复对
    byorg = collections.defaultdict(list)
    for s, r in smp.items():
        if r["SAMPLE_TYPE"] == "Metastasis" and r["METASTATIC_SITE"].strip():
            byorg[(pid(s), r["METASTATIC_SITE"].strip().lower())].append(r)
    noise = []
    for k, v in byorg.items():
        if len(v) >= 2:
            v.sort(key=lambda x: int(x["SAMPLE_COVERAGE"] or 0), reverse=True)
            noise.append((k[0], v[0], v[1]))

    cross_by = collections.defaultdict(list)
    for p, a, b in pairs:
        cross_by[a["CANCER_TYPE"]].append(int((a["sampleId"] in alt) != (b["sampleId"] in alt)))
    noise_by = collections.defaultdict(list)
    for p, a, b in noise:
        noise_by[smp[a["sampleId"]]["CANCER_TYPE"]].append(int((a["sampleId"] in alt) != (b["sampleId"] in alt)))

    common = [c for c in cross_by if len(noise_by.get(c, [])) >= 5]
    W = sum(len(cross_by[c]) for c in common)
    cross_w = sum(len(cross_by[c]) * np.mean(cross_by[c]) for c in common) / W
    noise_w = sum(len(cross_by[c]) * np.mean(noise_by[c]) for c in common) / W

    rng = np.random.default_rng(SEED)
    diffs = []
    for _ in range(2000):
        cc = nn = w = 0.0
        for c in common:
            ci = np.array(cross_by[c]); ni = np.array(noise_by[c])
            cc += len(ci) * ci[rng.integers(0, len(ci), len(ci))].mean()
            nn += len(ci) * ni[rng.integers(0, len(ni), len(ni))].mean()
            w += len(ci)
        diffs.append(100 * (cc - nn) / w)
    lo, hi = np.percentile(diffs, [2.5, 97.5])

    rows = [{
        "analysis": "tumour_type_standardised",
        "n_cancer_types": len(common),
        "pairs_covered": W,
        "pairs_total": len(pairs),
        "coverage_pct": round(100 * W / len(pairs), 1),
        "cross_pct_standardised": round(100 * cross_w, 2),
        "noise_pct_standardised": round(100 * noise_w, 2),
        "excess_pp": round(100 * (cross_w - noise_w), 2),
        "excess_lo": round(float(lo), 2),
        "excess_hi": round(float(hi), 2),
    }, {
        "analysis": "unstandardised_reference",
        "n_cancer_types": len(cross_by),
        "pairs_covered": len(pairs),
        "pairs_total": len(pairs),
        "coverage_pct": 100.0,
        "cross_pct_standardised": None,
        "noise_pct_standardised": None,
        "excess_pp": 7.9,
        "excess_lo": 4.8,
        "excess_hi": 11.0,
    }]
    out = os.path.join(TAB, "p1_standardised_excess.csv")
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w_ = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w_.writeheader(); w_.writerows(rows)
    print(f"cancer types standardised: {len(common)}; pairs covered {W}/{len(pairs)} ({100*W/len(pairs):.1f}%)")
    print(f"standardised cross {100*cross_w:.2f}%  noise {100*noise_w:.2f}%  excess {100*(cross_w-noise_w):.2f} pp (95% CI {lo:.2f} to {hi:.2f})")
    print("written:", out)


if __name__ == "__main__":
    main()
