#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
器官规则的把握度 / 精度边界分析（对应 analysis_plan_v2.0 的 P3/P4）

问题：以当前样本量，每个器官能"排除"多大的超额翻转？要给出可执行规则，还需要多少病例？

方法（全部为闭式解析，不模拟）：
  设 p0 = 该器官同器官重复翻转率（噪声地板），p1 = 原发→该器官翻转率，delta = p1 - p0
  SE(delta) = sqrt(p1(1-p1)/n1 + p0(1-p0)/n0)
  规则判定（v2.0 预设）：
    "可替代"      : 95% CI 上限 < 5 pp
    "不可替代"    : 95% CI 下限 > 5 pp
  由此反解所需样本量：
    可替代（真值 delta=0）      : 1.96*SE < 0.05  →  n
    不可替代（真值 delta=10pp） : delta - 1.96*SE > 0.05  →  n
输出：05_results/tables/p4_power_boundary.csv, 05_results/POWER_BOUNDARY_MEMO.md
"""
import os
import numpy as np
import pandas as pd

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "05_results")
TAB = os.path.join(OUT, "tables")
Z = 1.96
MARGIN = 0.05          # 规则判定的 5 个百分点边界
P0_STUDY = 0.1117      # 全研究噪声地板（62/555）

rule = pd.read_csv(os.path.join(TAB, "p4_organ_rule.csv"))

def se(p1, p0, n1, n0):
    return float(np.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0))

def n_for_substitutable(p0, ratio, margin=MARGIN):
    """真值 delta=0 时，要 CI 上限 < margin 所需 n1（n0 = ratio*n1）"""
    # Z*sqrt(p0(1-p0)/n1 + p0(1-p0)/(ratio*n1)) < margin
    return (Z / margin) ** 2 * (p0 * (1 - p0)) * (1 + 1 / ratio)

def n_for_not_substitutable(p0, ratio, delta=0.10, margin=MARGIN):
    """真值 delta=delta 时，要 CI 下限 > margin 所需 n1"""
    p1 = min(p0 + delta, 0.999)
    need_se = (delta - margin) / Z            # 要求 1.96*SE < delta-margin
    if need_se <= 0:
        return float("inf")
    return (Z / need_se) ** 2 * (p1 * (1 - p1) + p0 * (1 - p0) / ratio)

rows = []
for _, r in rule.iterrows():
    organ = r["organ"]
    n1 = int(r["cross_n"]); n0 = int(r["noise_n"]) if not pd.isna(r["noise_n"]) else 0
    p1 = r["cross_pct"] / 100
    p0_org = r["noise_pct"] / 100 if not pd.isna(r["noise_pct"]) else P0_STUDY
    if n0 < 15:                      # 噪声地板不足：无法估计该器官的重复变异
        rows.append({
            "organ": organ, "cross_n": n1, "noise_n": n0,
            "p_cross_pct": round(p1 * 100, 2), "p_noise_pct": round(p0_org * 100, 2),
            "current_CI_halfwidth_pp": None,
            "min_detectable_excess_pp": None,
            "can_exclude_excess_below_pp": None,
            "n_needed_substitutable_at_5pp": None,
            "n_needed_not_substitutable_at_10pp": None,
            "note": "噪声地板 n<15，无法评估",
        })
        continue
    s = se(p1, p0_org, n1, n0)
    hw = Z * s * 100                     # CI 半宽（pp）
    delta = (p1 - p0_org) * 100
    rows.append({
        "organ": organ, "cross_n": n1, "noise_n": n0,
        "p_cross_pct": round(p1 * 100, 2), "p_noise_pct": round(p0_org * 100, 2),
        "observed_excess_pp": round(delta, 2),
        "current_CI_halfwidth_pp": round(hw, 2),
        "min_detectable_excess_pp": round(hw, 2),
        "can_exclude_excess_below_pp": round(delta - hw, 2),
        "n_needed_substitutable_at_5pp": int(round(n_for_substitutable(p0_org, max(n0 / n1, 0.05)))),
        "n_needed_not_substitutable_at_10pp": int(round(n_for_not_substitutable(p0_org, max(n0 / n1, 0.05)))),
        "n_needed_if_noise_floor_known": int(round((Z / MARGIN) ** 2 * p0_org * (1 - p0_org))),
        "note": "",
    })

df = pd.DataFrame(rows).sort_values("cross_n", ascending=False)
df.to_csv(os.path.join(TAB, "p4_power_boundary.csv"), index=False, encoding="utf-8-sig")

# 全研究层面：当前能排除什么
n_all, k_all = 2278, int(round(0.191 * 2278))
se_all = se(435 / 2278, 62 / 555, 2278, 555)
lo_all = (435 / 2278 - 62 / 555 - Z * se_all) * 100
hi_all = (435 / 2278 - 62 / 555 + Z * se_all) * 100

lines = []
lines.append("# 器官规则的把握度边界（阶段 3 补充分析）\n")
lines.append(f"> 生成：{pd.Timestamp.now():%Y-%m-%d %H:%M}｜方法：解析法（闭式），非模拟\n")
lines.append("## 一、全研究层面：当前样本能排除什么\n")
lines.append(f"- 观测超额 **+7.9 pp**，95% CI **{lo_all:.1f}–{hi_all:.1f} pp**")
lines.append(f"- 因此可**排除**的区间是：器官归因超额 **< {max(lo_all,0):.1f} pp**")
lines.append(f"- 换言之：本研究足以证明「器官身份至少造成 {max(lo_all,0):.1f} 个百分点的额外翻转」，但**不足以证明「小于 5 pp 的器官效应不存在」**（那需要 CI 上限 <5 pp，当前上限 {hi_all:.1f} pp）\n")
lines.append("## 二、器官层面：每个器官现在能排除多少\n")
lines.append("| 器官 | 配对 n | 噪声 n | 观测超额 | 当前 CI 半宽 | 最小可检出超额 | 现在可排除的区间 | 判定为「可替代」(<5pp) 所需 n | 判定为「不可替代」(>5pp, 真值10pp) 所需 n |")
lines.append("|---|---|---|---|---|---|---|---|---|")
for _, r in df.iterrows():
    if r["current_CI_halfwidth_pp"] is None:
        lines.append(f"| {r['organ']} | {r['cross_n']} | {r['noise_n']} | — | — | — | — | — | — |")
        continue
    lines.append(
        f"| {r['organ']} | {r['cross_n']} | {r['noise_n']} | {r.get('observed_excess_pp')} pp | "
        f"±{r['current_CI_halfwidth_pp']} pp | {r['min_detectable_excess_pp']} pp | "
        f"{(('< ' + str(r['can_exclude_excess_below_pp']) + ' pp') if r['can_exclude_excess_below_pp'] > 0 else '无法排除（CI 下限 ≤0）')} | {r['n_needed_substitutable_at_5pp']} | "
        f"{r['n_needed_not_substitutable_at_10pp']} |")
lines.append("")
lines.append("## 三、结论（可直接写进 Discussion）\n")
feas = df[df["current_CI_halfwidth_pp"].notna()]
pos = feas[feas['can_exclude_excess_below_pp'] > 0]
lo_rng = f"{pos['can_exclude_excess_below_pp'].min():.1f}–{pos['can_exclude_excess_below_pp'].max():.1f}" if len(pos) else "无"
n_neg = len(feas) - len(pos)
n_sub = f"{int(feas['n_needed_substitutable_at_5pp'].min())}–{int(feas['n_needed_substitutable_at_5pp'].max())}"
n_not = f"{int(feas['n_needed_not_substitutable_at_10pp'].min())}–{int(feas['n_needed_not_substitutable_at_10pp'].max())}"
lines.append(f"1. **器官层面：只有 {len(pos)} 个器官（{', '.join(pos['organ']) if len(pos) else '无'}）的 CI 下限 >0，可排除「超额大于 {lo_rng} pp」；其余 {n_neg} 个器官的 CI 含 0，连「存在正超额」都无法确立**。全队列层面可排除「超额 <4.8 pp」。")
lines.append(f"2. 若要在单一器官上判定「可替代」（超额 <5 pp），按当前噪声地板与配对比例外推，约需 **{n_sub} 对**（现有多数器官为 190–419 对）。")
lines.append(f"3. 若要判定「不可替代」（超额 >5 pp，真值 10 pp），约需 **{n_not} 对**。")
bott = feas.loc[feas["n_needed_if_noise_floor_known"].idxmin()] if len(feas) else None
if bott is not None:
    lines.append(f"4. **瓶颈在噪声地板精度，不在跨界样本量**：若各器官的重复测量变异已被精确估计，判定「可替代」所需配对可降至约 {int(feas['n_needed_if_noise_floor_known'].min())}–{int(feas['n_needed_if_noise_floor_known'].max())} 对（现为 147–1710 对）。本数据集同器官重复仅 555 组（其中 503 组只有 2 个样本，全部组合也只有 671 对），**无法靠抽取方式补足**——必须引入第二个队列。")
lines.append("5. 因此本研究的正确表述是：**器官归因超额在全队列层面成立（+7.9 pp），但器官特异性规则超出本数据的精度边界**——这是精度限制，不是「没有差异」。")
with open(os.path.join(OUT, "POWER_BOUNDARY_MEMO.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")

print("\n".join(lines[:40]))
print("\n已写出 p4_power_boundary.csv 与 POWER_BOUNDARY_MEMO.md")
