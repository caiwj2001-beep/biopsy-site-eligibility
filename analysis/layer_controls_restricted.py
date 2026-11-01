#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""S3 补充表 Panel C：限制在 2,278 配对队列内的三层对照（敏感性）

与 layer_controls.py 使用同一套判定（30 基因、CNA>=2、同器官取覆盖度前二），
差别仅在于患者集合：这里只保留同时有原发与转移样本的患者。
输出：tables/layer_controls_restricted.csv
"""
import csv, collections, os

BASE = r"F:\数据库汇总\msk_impact_50k_2026"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "05_results", "tables")

E2S = {"1956":"EGFR","3845":"KRAS","673":"BRAF","2064":"ERBB2","4233":"MET","5290":"PIK3CA","672":"BRCA1",
       "675":"BRCA2","3417":"IDH1","3418":"IDH2","2260":"FGFR1","2263":"FGFR2","2261":"FGFR3","5727":"PTCH1",
       "6608":"SMO","207":"AKT1","2475":"MTOR","7248":"TSC1","7249":"TSC2","3815":"KIT","5156":"PDGFRA",
       "25":"ABL1","5979":"RET","238":"ALK","6098":"ROS1","4914":"NTRK1","2099":"ESR1","367":"AR",
       "1019":"CDK4","4193":"MDM2"}
GENES = set(E2S.values())
NONCODING = {"Silent","Intron","3'UTR","5'UTR","IGR","RNA","Flank","5'Flank","3'Flank"}

smp = {r["sampleId"]: r for r in csv.DictReader(open(os.path.join(BASE, "data_clinical_sample.csv"), encoding="utf-8-sig"))}
pid = lambda s: "-".join(s.split("-")[:2])
cov = lambda x: int(x["SAMPLE_COVERAGE"] or 0)

alt = collections.defaultdict(set)
for r in csv.DictReader(open(os.path.join(BASE, "data_mutations.csv"), encoding="utf-8-sig")):
    g = E2S.get(r["entrezGeneId"])
    if g and r["mutationType"] not in NONCODING:
        alt[r["sampleId"]].add(g)
with open(os.path.join(BASE, "data_cna.txt"), encoding="utf-8-sig") as f:
    head = f.readline().rstrip("\n").split("\t")[1:]
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] in GENES:
            for sid, v in zip(head, p[1:]):
                try:
                    if float(v) >= 2:
                        alt[sid].add(p[0])
                except ValueError:
                    pass
A = lambda s: s in alt

bypat = collections.defaultdict(list)
for s, r in smp.items():
    bypat[pid(s)].append(r)
index = {p for p, rs in bypat.items()
         if any(x["SAMPLE_TYPE"] == "Primary" for x in rs) and any(x["SAMPLE_TYPE"] == "Metastasis" for x in rs)}

def build(pats):
    same, cross = [], []
    for p in pats:
        d = collections.defaultdict(list)
        for x in bypat[p]:
            if x["SAMPLE_TYPE"] == "Metastasis" and x["METASTATIC_SITE"].strip():
                d[x["METASTATIC_SITE"].strip().lower()].append(x)
        for o, v in d.items():
            if len(v) >= 2:
                v = sorted(v, key=cov, reverse=True)
                same.append((p, v[0], v[1]))
        if len(d) >= 2:
            o = sorted(d)
            a = sorted(d[o[0]], key=cov, reverse=True)[0]
            b = sorted(d[o[1]], key=cov, reverse=True)[0]
            cross.append((p, a, b))
    return same, cross

rows = []
for label, pats in [("all_patients", set(bypat)), ("restricted_to_paired_cohort", index)]:
    s_, c_ = build(pats)
    ds = sum(1 for _, a, b in s_ if A(a["sampleId"]) != A(b["sampleId"]))
    dc = sum(1 for _, a, b in c_ if A(a["sampleId"]) != A(b["sampleId"]))
    rows.append({"cohort": label,
                 "same_organ_pairs": len(s_), "same_organ_discordant": ds,
                 "same_organ_pct": round(100 * ds / len(s_), 2) if s_ else "",
                 "cross_organ_pairs": len(c_), "cross_organ_discordant": dc,
                 "cross_organ_pct": round(100 * dc / len(c_), 2) if c_ else ""})

out = os.path.join(TAB, "layer_controls_restricted.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
print("written:", out)
