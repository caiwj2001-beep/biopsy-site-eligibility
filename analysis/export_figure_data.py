#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""导出绘图所需的患者级数据（OS KM 与器官级），供 04_code/figures.R 使用。"""
import csv, os
from collections import defaultdict

BASE = r"F:\数据库汇总\msk_impact_50k_2026"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "05_results", "tables")
os.makedirs(OUT, exist_ok=True)

TARGET = {"1956","3845","673","2064","4233","5290","672","675","3417","3418","2260","2263","2261",
          "5727","6608","207","2475","7248","7249","3815","5156","25","5979","238","6098","4914",
          "2099","367","1019","4193"}
NONCODING = {"Silent","Intron","3'UTR","5'UTR","IGR","RNA","Flank","5'Flank","3'Flank"}

samples = {r["sampleId"]: r for r in csv.DictReader(open(os.path.join(BASE,"data_clinical_sample.csv"), encoding="utf-8-sig"))}
patients = {}
hdr = None
for line in open(os.path.join(BASE,"data_clinical_patient.txt"), encoding="utf-8-sig"):
    if line.startswith("#"): continue
    p = line.rstrip("\n").split("\t")
    if p[0] == "PATIENT_ID": hdr = p; continue
    patients[p[0]] = dict(zip(hdr, p))

pid = lambda s: "-".join(s.split("-")[:2])
alt = defaultdict(set)
for r in csv.DictReader(open(os.path.join(BASE,"data_mutations.csv"), encoding="utf-8-sig")):
    if r["entrezGeneId"] in TARGET and r["mutationType"] not in NONCODING:
        alt[r["sampleId"]].add(r["entrezGeneId"])
for line in open(os.path.join(BASE,"data_cna.txt"), encoding="utf-8-sig"):
    parts = line.rstrip("\n").split("\t")
    if parts[0] not in {"EGFR","ERBB2","MET","MDM2","CDK4","FGFR1","FGFR2","FGFR3","KIT","PIK3CA"}: continue
    for sid, v in zip(parts[1:], parts[1:]):
        pass
# CNA 已并入主分析；此处 OS 图只需 discordant 标签，用与主分析一致的合并定义
cna_rows = {}
with open(os.path.join(BASE,"data_cna.txt"), encoding="utf-8-sig") as f:
    head = f.readline().rstrip("\n").split("\t")[1:]
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] not in {"EGFR","ERBB2","MET","MDM2","CDK4","FGFR1","FGFR2","FGFR3","KIT","PIK3CA","AR","ESR1"}: continue
        for sid, v in zip(head, p[1:]):
            try:
                if float(v) >= 2: alt[sid].add(p[0])
            except ValueError:
                pass

bypat = defaultdict(list)
for s, r in samples.items(): bypat[pid(s)].append(r)
rows = []
for p, rs in bypat.items():
    prim = [x for x in rs if x["SAMPLE_TYPE"]=="Primary"]
    met  = [x for x in rs if x["SAMPLE_TYPE"]=="Metastasis"]
    if not (prim and met): continue
    prim.sort(key=lambda x:int(x["SAMPLE_COVERAGE"] or 0), reverse=True)
    met.sort(key=lambda x:int(x["SAMPLE_COVERAGE"] or 0), reverse=True)
    pr, mt = prim[0], met[0]
    pt = patients.get(p, {})
    try: t = float(pt["OS_MONTHS"])
    except (KeyError, ValueError): continue
    ev = 1 if pt.get("OS_STATUS","").strip().upper()=="DECEASED" else 0
    disc = int((pr["sampleId"] in alt) != (mt["sampleId"] in alt))
    rows.append({"patient":p, "os_months":t, "event":ev, "discordant":disc,
                 "cancer":pr["CANCER_TYPE"], "organ":(mt["METASTATIC_SITE"] or "Unspecified")})

with open(os.path.join(OUT,"os_km_input.csv"),"w",newline="",encoding="utf-8-sig") as f:
    w=csv.DictWriter(f, fieldnames=["patient","os_months","event","discordant","cancer","organ"])
    w.writeheader(); w.writerows(rows)
print("os_km_input rows:", len(rows), "events:", sum(r['event'] for r in rows))
