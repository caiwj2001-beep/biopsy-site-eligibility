#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
第二队列（独立）复现：AURORA US Network 转移性乳腺癌（cBioPortal brca_aurora_2023）

背景：本项目的核心终点需要"同一患者的原发灶 + 转移灶"配对。公开数据中，
      cBioPortal 551 个研究中唯一非 MSK 的配对队列即 AURORA（48 对）。
      本脚本用与主分析完全相同的 30 基因可靶向定义，计算该独立队列的翻转率，
      并与 MSK 乳腺癌亚组（322 对）比较。

注意（必须写进论文）：
  - AURORA 为 WES，无 CNA 扩增判定（本研究只用 SNV/indel，与主分析的"完整定义"不完全同构）；
  - 仅乳腺癌一个瘤种，n=48，只能做方向一致性判断，不能做等效性检验；
  - 该队列与本项目主队列无样本重叠（不同机构、不同测序平台）。
"""
import json, os, collections, urllib.request
import pandas as pd

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "05_results")
TAB = os.path.join(OUT, "tables"); os.makedirs(TAB, exist_ok=True)
SID = "brca_aurora_2023"
API = "https://www.cbioportal.org/api"
TARGET = ["EGFR","KRAS","BRAF","ERBB2","MET","PIK3CA","BRCA1","BRCA2","IDH1","IDH2",
          "FGFR1","FGFR2","FGFR3","PTCH1","SMO","AKT1","MTOR","TSC1","TSC2","KIT",
          "PDGFRA","ABL1","RET","ALK","ROS1","NTRK1","ESR1","AR","CDK4","MDM2"]
NONCODING = {"Silent","Intron","3'UTR","5'UTR","IGR","RNA","Flank","5'Flank","3'Flank"}

def get(u):
    return json.load(urllib.request.urlopen(u, timeout=120))

def post(u, payload):
    req = urllib.request.Request(u, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=300))

def wilson(k, n, z=1.96):
    if n == 0: return (float("nan"),)*3
    p = k/n; d = 1+z*z/n; c = p+z*z/(2*n)
    m = z*((p*(1-p)/n + z*z/(4*n*n))**0.5)
    return (p*100, (c-m)/d*100, (c+m)/d*100)

def main():
    clin = get(f"{API}/studies/{SID}/clinical-data?clinicalDataType=SAMPLE&attributeId=SAMPLE_TYPE&projection=SUMMARY")
    sid2pat = {x["sampleId"]: x["patientId"] for x in clin}
    sid2type = {x["sampleId"]: (x["value"] or "").lower() for x in clin}
    print(f"samples={len(sid2type)} patients={len(set(sid2pat.values()))}")

    pairs = []
    bypat = collections.defaultdict(dict)
    for s, p in sid2pat.items():
        t = sid2type[s]
        if "primary" in t: bypat[p].setdefault("primary", []).append(s)
        elif "metasta" in t: bypat[p].setdefault("met", []).append(s)
    for p, d in bypat.items():
        if d.get("primary") and d.get("met"):
            pairs.append((p, sorted(d["primary"])[0], sorted(d["met"])[0]))
    print(f"paired patients (primary + metastasis) = {len(pairs)}")

    samples = sorted({s for _, a, b in pairs for s in (a, b)})
    alt = collections.defaultdict(set)
    for i in range(0, len(samples), 40):
        batch = samples[i:i+40]
        recs = post(f"{API}/molecular-profiles/{SID}_mutations/mutations/fetch", {"sampleIds": batch})
        for r in recs:
            g = r.get("hugoGeneSymbol") or r.get("keyword", "").split(" ")[0]
            if g in TARGET and r.get("mutationType") not in NONCODING:
                alt[r["sampleId"]].add(g)

    flags = []
    for p, a, b in pairs:
        flags.append((a in alt) != (b in alt))
    k = sum(flags); n = len(flags)
    rate = wilson(k, n)

    rows = [{
        "cohort": "AURORA_US_breast_WES", "platform": "WES (SNV/indel only)",
        "pairs": n, "discordant": k, "discordance_pct": round(rate[0], 2),
        "ci_lo": round(rate[1], 2), "ci_hi": round(rate[2], 2),
        "note": "independent of MSK; single tumour type; no CNA amplification",
    }]
    pd.DataFrame(rows).to_csv(os.path.join(TAB, "replication_aurora.csv"), index=False, encoding="utf-8-sig")
    print(f"AURORA: {k}/{n} = {rate[0]:.1f}% (95%CI {rate[1]:.1f}-{rate[2]:.1f})")

main()
