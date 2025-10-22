#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combine attention scores from multiple runs into a single ranked table,
then compute per-edge confidence levels based on cross-celltype rank consistency.
Example:
python combine_attention_scores.py \
  --epoch 1000 \
  --root-dir ./experiments/multi_run \
  --celltype-file ./symbol_pair.txt \
  --output ./attention_layer1_epoch1000_combined_rank.csv \
  --output-wide ./final_attention_layer1_result.csv
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Combine attention scores across multiple runs.")
    parser.add_argument("--root-dir", type=str, default="./experiments/multi_run",
                        help="Directory that contains run_XX subdirectories.")
    parser.add_argument("--n-runs", type=int, default=10,
                        help="Number of runs to combine (expects run_01 … run_N).")
    parser.add_argument("--celltype-file", type=str, default="./symbol_pairs.txt",
                        help="Tab-delimited file with cell-type info; must contain 'celltype_stim' column.")
    parser.add_argument("--epoch", type=int, required=True,
                        help="Epoch number – file name is attention_layer1_epoch{epoch}.csv.")
    parser.add_argument("--output", type=str, default="./all_attention_layer1_combined_rank.csv",
                        help="Output path for the combined CSV (directories created if needed).")
    parser.add_argument("--output-wide", type=str, default="./all_attention_layer1_wide_confidence.csv",
                        help="Output path for the wide-format confidence table.")
    return parser.parse_args()


def load_run_scores(root: Path, run_name: str, csv_name: str) -> pd.DataFrame:
    csv_path = root / run_name / csv_name
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing file: {csv_path}")

    df_raw = pd.read_csv(csv_path, header=0, usecols=[0])
    series = pd.to_numeric(df_raw.iloc[:, 0], errors="coerce")
    series.name = run_name
    return series.to_frame()


def main() -> None:
    args = parse_args()
    csv_name = f"attention_layer1_epoch{args.epoch}.csv"

    root_dir = Path(args.root_dir).expanduser().resolve()
    if not root_dir.is_dir():
        sys.exit(f"[ERROR] Root directory not found: {root_dir}")

    # === 1️加载多个 run 的注意力得分 ===
    run_names = [f"run_{i:02d}" for i in range(1, args.n_runs + 1)]
    dfs = []
    for run in run_names:
        try:
            dfs.append(load_run_scores(root_dir, run, csv_name))
        except FileNotFoundError as exc:
            sys.exit(str(exc))

    df_all = pd.concat(dfs, axis=1)

    # === 2️ 加载 celltype 注释文件 ===
    celltype_path = Path(args.celltype_file).expanduser().resolve()
    if not celltype_path.exists():
        sys.exit(f"[ERROR] Cell-type file not found: {celltype_path}")
    celltype_df = pd.read_csv(celltype_path, sep="\t", header=0)

    # === 3️ 合并信息 ===
    df = pd.concat([celltype_df, df_all], axis=1)
    score_cols = run_names
    n_runs = len(score_cols)
    group_key = "celltype_stim"

    # === 4️各种组内排名指标 ===
    df["mean_score"] = df[score_cols].mean(axis=1, skipna=True)
    df["rank_mean_in_grp"] = df.groupby(group_key)["mean_score"].rank(method="min", ascending=False)

    rank_df = df.groupby(group_key)[score_cols].rank(method="min", ascending=False)
    df["rank_product"] = np.exp(np.log(rank_df).sum(axis=1) / n_runs)
    df["rank_rp_in_grp"] = df.groupby(group_key)["rank_product"].rank(method="min", ascending=True)

    eps = np.finfo(float).eps
    geo_mean = np.exp(np.log(df[score_cols].clip(lower=eps)).sum(axis=1) / n_runs)
    df["geo_mean_score"] = geo_mean
    df["rank_gm_in_grp"] = df.groupby(group_key)["geo_mean_score"].rank(method="min", ascending=False)

    # === 5️输出初步的合并文件 ===
    output_path = Path(args.output).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"[INFO] Combined table written to {output_path}")

    # === 6️R 部分逻辑的 Python 实现：宽表 + rank consistency ===
    print("[INFO] Computing per-edge cross-celltype rank consistency...")

    df_wide = df.pivot_table(index="edge", columns="celltype_stim", values="rank_rp_in_grp")
    df_wide = df_wide.reset_index()

    # rowMins → 每行最小 rank
    df_wide["rank"] = df_wide.iloc[:, 1:].min(axis=1, skipna=True)
    df_wide["final_rank"] = df_wide["rank"].rank(method="min", ascending=True)

    # 计算 80% 分位点阈值
    threshold = df_wide["final_rank"].quantile(0.8)

    # 新增置信度列
    df_wide["confidence"] = np.where(
        df_wide["final_rank"] <= threshold, "High_confidence", "Low_confidence"
    )

    # 输出最终宽表
    output_wide_path = Path(args.output_wide).expanduser().resolve()
    output_wide_path.parent.mkdir(parents=True, exist_ok=True)
    df_wide.to_csv(output_wide_path, index=False)
    print(f"[INFO] Wide-format confidence table written to {output_wide_path}")


if __name__ == "__main__":
    main()
