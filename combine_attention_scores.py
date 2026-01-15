#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch_geometric.data import Data

from GATv2Conv_CellNEST import GATv2Conv

warnings.filterwarnings("ignore")


# =====================================================
# Encoder
# =====================================================
class Encoder(nn.Module):
    def __init__(self, in_channels, hidden_channels, edge_dim, heads=1):
        super().__init__()

        self.gat1 = GATv2Conv(
            in_channels,
            hidden_channels,
            edge_dim=edge_dim,
            heads=heads,
            concat=False,
            add_self_loops=False,
        )
        self.bn1 = nn.BatchNorm1d(hidden_channels)
        self.dropout1 = nn.Dropout(0.3)

        self.gat2 = GATv2Conv(
            hidden_channels,
            hidden_channels,
            edge_dim=edge_dim,
            heads=heads,
            concat=False,
            add_self_loops=False,
        )
        self.bn2 = nn.BatchNorm1d(hidden_channels)
        self.dropout2 = nn.Dropout(0.3)

        self.prelu = nn.PReLU(hidden_channels)
        self.attn1_raw = None

    def forward(self, x, edge_index, edge_attr=None):
        x, _, attn_raw = self.gat1(
            x,
            edge_index,
            edge_attr=edge_attr,
            return_attention_weights=True,
        )
        self.attn1_raw = attn_raw

        x = self.dropout1(self.bn1(x))
        x, _, _ = self.gat2(
            x,
            edge_index,
            edge_attr=edge_attr,
            return_attention_weights=True,
        )
        x = self.dropout2(self.bn2(x))
        return self.prelu(x)


# =====================================================
# Args
# =====================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract GATv2 attention scores from multiple runs and aggregate them."
    )

    parser.add_argument("--node-emb", required=True, help="Node embedding file")
    parser.add_argument("--edge-pairs", required=True, help="Edge index pairs file")
    parser.add_argument("--edge-feats", required=True, help="Edge feature file")
    parser.add_argument("--symbol-pairs", required=True, help="Symbol/celltype metadata file")

    parser.add_argument(
        "--weight-dir",
        required=True,
        help="Directory containing run_XX subfolders (e.g. ./experiments/multi_run)",
    )

    parser.add_argument("--n-runs", type=int, default=10, help="Number of runs")
    parser.add_argument("--in-dim", type=int, default=512, help="Input node dimension")
    parser.add_argument("--hidden-dim", type=int, default=64, help="Hidden dimension")
    parser.add_argument("--heads", type=int, default=1, help="GAT heads")

    parser.add_argument(
        "--group-key",
        default="celltype_stim",
        help="Grouping column in metadata",
    )

    parser.add_argument(
        "--output",
        default="all_attention_combined_with_rank.csv",
        help="Output CSV file",
    )

    return parser.parse_args()


# =====================================================
# Main
# =====================================================
def main():
    args = parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[INFO] Device: {device}")

    # -----------------------------
    # Load graph data
    # -----------------------------
    node_emb = torch.tensor(
        pd.read_csv(args.node_emb, sep="\t", header=None).values,
        dtype=torch.float,
    )

    edge_pairs = torch.tensor(
        pd.read_csv(args.edge_pairs, sep="\t", header=None).values,
        dtype=torch.long,
    )

    edge_feats = torch.tensor(
        pd.read_csv(args.edge_feats, sep="\t", header=None).values,
        dtype=torch.float,
    )

    edge_dim = edge_feats.shape[1]

    data = Data(
        x=node_emb,
        edge_index=edge_pairs.t(),
        edge_attr=edge_feats,
    ).to(device)

    # -----------------------------
    # Extract attention scores
    # -----------------------------
    attn_scores_df = pd.DataFrame()

    for run_id in range(1, args.n_runs + 1):
        run_dir = Path(args.weight_dir) / f"run_{run_id:02d}"
        weight_path = run_dir / "encoder_epoch1000.pt"

        if not weight_path.exists():
            raise FileNotFoundError(f"Missing weight file: {weight_path}")

        print(f"[INFO] Run {run_id}: loading {weight_path}")

        model = Encoder(
            in_channels=args.in_dim,
            hidden_channels=args.hidden_dim,
            edge_dim=edge_dim,
            heads=args.heads,
        ).to(device)

        model.load_state_dict(torch.load(weight_path, map_location=device))
        model.eval()

        with torch.no_grad():
            model(data.x, data.edge_index, data.edge_attr)

        attn_score = (
            torch.sigmoid(model.attn1_raw)
            .detach()
            .cpu()
            .squeeze(-1)
            .numpy()
        )

        attn_scores_df[f"attn_score{run_id}"] = attn_score

        del model
        torch.cuda.empty_cache()

    # -----------------------------
    # Merge metadata
    # -----------------------------
    meta_df = pd.read_csv(args.symbol_pairs, sep="\t", header=0)
    df = pd.concat([meta_df, attn_scores_df], axis=1)

    score_cols = attn_scores_df.columns.tolist()
    n_runs = len(score_cols)
    group_key = args.group_key

    # -----------------------------
    # Ranking
    # -----------------------------
    df["mean_score"] = df[score_cols].mean(axis=1, skipna=True)
    df["rank_mean_in_grp"] = (
        df.groupby(group_key)["mean_score"]
        .rank(method="min", ascending=False)
    )

    rank_df = df.groupby(group_key)[score_cols].rank(
        method="min", ascending=False
    )
    df["rank_product"] = np.exp(
        np.log(rank_df).sum(axis=1) / n_runs
    )
    df["rank_rp_in_grp"] = (
        df.groupby(group_key)["rank_product"]
        .rank(method="min", ascending=True)
    )

    eps = np.finfo(float).eps
    df["geo_mean_score"] = np.exp(
        np.log(df[score_cols].clip(lower=eps)).sum(axis=1) / n_runs
    )
    df["rank_gm_in_grp"] = (
        df.groupby(group_key)["geo_mean_score"]
        .rank(method="min", ascending=False)
    )

    # -----------------------------
    # Save
    # -----------------------------
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    print(f"[INFO] Saved result to: {output_path}")


if __name__ == "__main__":
    main()

