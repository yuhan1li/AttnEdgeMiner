#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Train a two‑layer GATv2 encoder + DGI model for link‑prediction‑style
attention mining **multiple times** and aggregate the results.

Example (3 runs):
    python train_gatv_2_dgi_multi.py \
        --edge-pairs-true ./index_pairs_true.txt \
        --edge-feats-true ./edge_feature_true.txt \
        --edge-pairs ./index_pairs_all.txt \
        --edge-feats ./edge_feature_all.txt \
        --node-emb  ./gene_embedding.txt \
        --pos-edges ./Positive_edge.txt \
        --neg-edges ./Negative_edge.txt \
        --outdir     ./experiments/multi_run \
        --epochs     1000 \
        --runs       10 \
        --dgi-scale  0.1

Outputs
-------
<outdir>/
    run_01/  …     # artifacts from each individual run
    run_02/  …
    run_03/  …
    all_runs_history.csv         # concatenated history of all runs
    mean_history.csv             # epoch‑wise mean curves
    mean_curves.png              # mean loss / metric curves
    final_metrics_summary.csv    # table of final‑epoch metrics for every run
"""
from __future__ import annotations
import argparse
from pathlib import Path
import warnings
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from torch_geometric.data import Data
from torch_geometric.nn import DeepGraphInfomax

# -----------------------------------------------------------------------------
# 0) Suppress deprecation warnings
# -----------------------------------------------------------------------------
warnings.filterwarnings("ignore", category=UserWarning, message=r"nn\.functional\.tanh is deprecated")

# -----------------------------------------------------------------------------
# 1) Import user‑provided GATv2Conv implementation
# -----------------------------------------------------------------------------
try:
    from GATv2Conv_CellNEST import GATv2Conv  # pyright: ignore
except ImportError as err:
    raise SystemExit("Cannot import GATv2Conv from GATv2Conv_CellNEST: " + str(err))

# -----------------------------------------------------------------------------
# 2) Loss function: binary focal loss on logits
# -----------------------------------------------------------------------------

def focal_loss_with_logits(
    logits: torch.Tensor,
    targets: torch.Tensor,
    *,
    alpha: float = 0.25,
    gamma: float = 2.0,
    reduction: str = "mean",
):
    """Binary focal loss built on BCEWithLogits."""
    bce = F.binary_cross_entropy_with_logits(logits, targets.float(), reduction="none")
    pt = torch.exp(-bce)
    focal = (1.0 - pt) ** gamma
    loss = alpha * focal * bce
    if reduction == "mean":
        return loss.mean()
    if reduction == "sum":
        return loss.sum()
    return loss

# -----------------------------------------------------------------------------
# 3) Encoder + DGI wrapper
# -----------------------------------------------------------------------------

class Encoder(nn.Module):
    """Two-layer GATv2 encoder that stores raw attention logits."""

    def __init__(self, in_channels: int, hidden_channels: int, edge_dim: int, heads: int = 1):
        super().__init__()
        self.edge_dim = edge_dim  # 自动识别后传入

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
        self.attn1_raw: torch.Tensor | None = None  # 存储 raw attention logits

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor, edge_attr: torch.Tensor | None = None):
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


class DGIWrapper(nn.Module):
    """Encoder + Deep Graph Infomax."""

    def __init__(self, in_channels: int, hidden_channels: int, edge_dim: int, heads: int, dropout: float):
        super().__init__()
        self.encoder = Encoder(in_channels, hidden_channels, edge_dim, heads)
        self.dgi = DeepGraphInfomax(
            hidden_channels,
            self.encoder,
            corruption=self._corrupt,
            summary=lambda z, *_: torch.sigmoid(z.mean(dim=0)),
        )

    def _corrupt(self, x: torch.Tensor, edge_index: torch.Tensor, edge_attr: torch.Tensor | None):
        return x[torch.randperm(x.size(0))], edge_index, edge_attr

    def forward(self, data: Data):
        return self.dgi(data.x, data.edge_index, data.edge_attr)

# -----------------------------------------------------------------------------
# 4) Utilities
# -----------------------------------------------------------------------------

def sample_edge_indices(pos_edges: torch.Tensor, neg_edges: torch.Tensor, *, num: int, device: torch.device):
    """Return indices of *num* positive & negative edges and their labels."""
    num = min(num, pos_edges.size(0), neg_edges.size(0))
    pos_idx = pos_edges[torch.randperm(pos_edges.size(0))[:num], 2]
    neg_idx = neg_edges[torch.randperm(neg_edges.size(0))[:num], 2]
    combined = torch.cat([pos_idx, neg_idx], dim=0).to(device)
    labels = torch.cat([torch.ones(num), torch.zeros(num)], dim=0).to(device)
    return combined, labels

# -----------------------------------------------------------------------------
# 5) Single‑run training function (was `train`)
# -----------------------------------------------------------------------------

def train_single_run(
    *,
    edge_pairs_true_path: Path,
    edge_feats_true_path: Path,
    edge_pairs_path: Path,
    edge_feats_path: Path,
    node_emb_path: Path,
    pos_edges_path: Path,
    neg_edges_path: Path,
    outdir: Path,
    epochs: int = 1000,
    hidden_dim: int = 64,
    heads: int = 1,
    sample_size: int = 2000,
    lr: float = 1e-3,
    dgi_scale: float = 0.1,
) -> pd.DataFrame:
    """Train once and return the history DataFrame."""

    outdir.mkdir(parents=True, exist_ok=True)

    # ---------------- 1) Data ----------------
    edge_pairs_true = torch.tensor(pd.read_csv(edge_pairs_true_path, sep="\t", header=None).values, dtype=torch.long)
    edge_feats_true = torch.tensor(pd.read_csv(edge_feats_true_path, sep="\t", header=None).values, dtype=torch.float)
    edge_pairs = torch.tensor(pd.read_csv(edge_pairs_path, sep="\t", header=None).values, dtype=torch.long)
    edge_feats = torch.tensor(pd.read_csv(edge_feats_path, sep="\t", header=None).values, dtype=torch.float)
    node_emb = torch.tensor(pd.read_csv(node_emb_path, sep="\t", header=None).values, dtype=torch.float)
    pos_edges = torch.tensor(pd.read_csv(pos_edges_path, sep="\t", header=None).values, dtype=torch.long)
    neg_edges = torch.tensor(pd.read_csv(neg_edges_path, sep="\t", header=None).values, dtype=torch.long)
    edge_dim = edge_feats.shape[1]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data_real = Data(x=node_emb, edge_index=edge_pairs_true.t(), edge_attr=edge_feats_true).to(device)
    data_with_fake = Data(x=node_emb, edge_index = edge_pairs.t(), edge_attr = edge_feats).to(device)
    
    # ---------------- 2) Model ----------------
    model = DGIWrapper(node_emb.size(1), hidden_dim, edge_dim, heads, dropout=0.2).to(device)
    enc = model.encoder
    optim = torch.optim.Adam(model.parameters(), lr=lr)

    # ---------------- 3) History ----------------
    hist: dict[str, list[float]] = {
        "loss1": [],
        "loss_dgi": [],
        "loss_total": [],
        "accuracy": [],
        "precision": [],
        "recall": [],
        "f1": [],
    }

    for epoch in range(1, epochs + 1):
        enc.train(); model.train()

        # Forward --------------------------
        _ = enc(data_with_fake.x, data_with_fake.edge_index, data_with_fake.edge_attr)
        logits_raw = enc.attn1_raw.mean(-1)  # [E]

        edge_sel, labels = sample_edge_indices(pos_edges, neg_edges, num=sample_size, device=device)
        logits_sel = logits_raw[edge_sel]

        loss1 = focal_loss_with_logits(logits_sel, labels, alpha=0.75, gamma=1.5)
        pos_z, neg_z, summary = model(data_real)
        loss_dgi = model.dgi.loss(pos_z, neg_z, summary)
        total_loss = loss1 + dgi_scale * loss_dgi

        # Backprop ------------------------
        optim.zero_grad(); total_loss.backward(); optim.step()

        # Metrics --------------------------
        probs = torch.sigmoid(logits_sel).detach().cpu().numpy()
        preds = (probs >= 0.5).astype(float)
        labels_np = labels.cpu().numpy()
        acc = accuracy_score(labels_np, preds)
        prec = precision_score(labels_np, preds)
        rec = recall_score(labels_np, preds)
        f1 = f1_score(labels_np, preds)

        # Log ------------------------------
        hist["loss1"].append(loss1.item())
        hist["loss_dgi"].append(dgi_scale * loss_dgi.item())
        hist["loss_total"].append(total_loss.item())
        hist["accuracy"].append(acc)
        hist["precision"].append(prec)
        hist["recall"].append(rec)
        hist["f1"].append(f1)

        if epoch % 20 == 0 or epoch == 1:
            print(
                f"Ep {epoch:04d} | total {total_loss:.4f} | link {loss1:.4f} | dgi {loss_dgi:.4f} "
                f"| Acc {acc:.4f} | Prec {prec:.4f} | Rec {rec:.4f} | F1 {f1:.4f}")

    # ---------------- 4) Save per‑run artifacts ----------------
    df_hist = pd.DataFrame({"epoch": range(1, epochs + 1), **hist})
    df_hist.to_csv(outdir / "training_history.csv", index=False)
    print(f"[Saved] history -> {outdir/'training_history.csv'}")

    # Curves ----------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(df_hist.epoch, df_hist.loss1, label="loss1")
    ax1.plot(df_hist.epoch, df_hist.loss_dgi, label="loss_dgi")
    ax1.plot(df_hist.epoch, df_hist.loss_total, label="total", linewidth=2)
    ax1.set_xlabel("Epoch"); ax1.set_ylabel("Loss"); ax1.set_title("Loss Curves"); ax1.legend()

    ax2.plot(df_hist.epoch, df_hist.accuracy, label="Acc")
    ax2.plot(df_hist.epoch, df_hist.precision, label="Prec")
    ax2.plot(df_hist.epoch, df_hist.recall, label="Rec")
    ax2.plot(df_hist.epoch, df_hist.f1, label="F1", linewidth=2)
    ax2.set_xlabel("Epoch"); ax2.set_ylabel("Score"); ax2.set_title("Metrics Curves"); ax2.legend()

    plt.tight_layout(); plt.savefig(outdir / "loss_metrics.png", dpi=300); plt.close()

    # Attention scores -----------------------------------------
    raw_attn = enc.attn1_raw.detach().cpu().squeeze(-1).numpy()
    attn = torch.sigmoid(enc.attn1_raw).detach().cpu().squeeze(-1).numpy()
    # 保存未 sigmoid 的注意力分数
    pd.DataFrame({"attention": raw_attn}).to_csv(outdir / f"attention_layer1_epoch{epochs}_raw.csv", index=False)
    pd.DataFrame({"attention": attn}).to_csv(outdir / f"attention_layer1_epoch{epochs}.csv", index=False)
    
    torch.save(enc.state_dict(), outdir / f"encoder_epoch{epochs}.pt")
    print(f"[Saved] model -> {outdir / f'encoder_epoch{epochs}.pt'}")

    return df_hist

# -----------------------------------------------------------------------------
# 6) Aggregation helpers
# -----------------------------------------------------------------------------

def aggregate_runs(base_outdir: Path, runs: int, epochs: int):
    """Read perrun histories, concatenate, and save mean curves & summary."""
    dfs = []
    for i in range(1, runs + 1):
        csv = base_outdir / f"run_{i:02d}" / "training_history.csv"
        if not csv.exists():
            print(f"[WARN] Missing history: {csv}")
            continue
        df = pd.read_csv(csv)
        df["run"] = i
        dfs.append(df)
    if not dfs:
        print("[WARN] No histories to aggregate."); return

    all_hist = pd.concat(dfs, ignore_index=True)
    all_hist.to_csv(base_outdir / "all_runs_history.csv", index=False)

    # Epoch‑wise mean curves
    mean_hist = all_hist.groupby("epoch").mean(numeric_only=True).reset_index()
    mean_hist.to_csv(base_outdir / "mean_history.csv", index=False)

    # Plot mean curves
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    ax1.plot(mean_hist.epoch, mean_hist.loss1, label="loss1")
    ax1.plot(mean_hist.epoch, mean_hist.loss_dgi, label="loss_dgi")
    ax1.plot(mean_hist.epoch, mean_hist.loss_total, label="total", linewidth=2)
    ax1.set_title("Mean Loss Curves"); ax1.set_xlabel("Epoch"); ax1.set_ylabel("Loss"); ax1.legend()

    ax2.plot(mean_hist.epoch, mean_hist.accuracy, label="Acc")
    ax2.plot(mean_hist.epoch, mean_hist.precision, label="Prec")
    ax2.plot(mean_hist.epoch, mean_hist.recall, label="Rec")
    ax2.plot(mean_hist.epoch, mean_hist.f1, label="F1", linewidth=2)
    ax2.set_title("Mean Metrics"); ax2.set_xlabel("Epoch"); ax2.set_ylabel("Score"); ax2.legend()

    plt.tight_layout(); plt.savefig(base_outdir / "mean_curves.png", dpi=300); plt.close()

    # Final‑epoch summary table --------------------------------
    final_rows = all_hist[all_hist.epoch == epochs]
    final_rows.to_csv(base_outdir / "final_metrics_summary.csv", index=False)
    print(f"[Saved] Aggregated artifacts in {base_outdir}")

# -----------------------------------------------------------------------------
# 7) CLI
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Train a two-layer GATv2 + DGI encoder for link-prediction attention mining "
            "multiple times and aggregate the results."
        )
    )
    p.add_argument(
        "--edge-pairs-true", required=True, type=Path,
        help="TSV file of REAL graph edges (src, dst, optional idx). Used for DGI training."
    )
    p.add_argument(
        "--edge-feats-true", required=True, type=Path,
        help="Edge feature matrix (same row order as --edge-pairs-true). Used as edge_attr in real graph."
    )
    p.add_argument(
        "--edge-pairs", required=True, type=Path,
        help="TSV of ALL edges (true + fake) to compute attention scores on."
    )
    p.add_argument(
        "--edge-feats", required=True, type=Path,
        help="Edge features for --edge-pairs. Shape: [num_edges, edge_dim]."
    )
    p.add_argument(
        "--node-emb", required=True, type=Path,
        help="Initial node embeddings (gene/cell vectors). Shape: [num_nodes, dim]."
    )
    p.add_argument(
        "--pos-edges", required=True, type=Path,
        help="Known positive edges (src, dst, edge_index) for supervised link loss."
    )
    p.add_argument(
        "--neg-edges", required=True, type=Path,
        help="Negative (non-existing) edges used as contrastive samples."
    )
    p.add_argument(
        "--outdir", required=True, type=Path,
        help="Base directory to save outputs of all runs."
    )
    p.add_argument("--epochs", type=int, default=1000,
                   help="Training epochs for each run (default = 1000).")
    p.add_argument("--hidden-dim", type=int, default=64,
                   help="Hidden dimension of GATv2 layers.")
    p.add_argument("--heads", type=int, default=1,
                   help="Number of attention heads.")
    p.add_argument("--sample-size", type=int, default=2000,
                   help="Number of pos/neg edges sampled per epoch for link loss.")
    p.add_argument("--lr", type=float, default=1e-3,
                   help="Learning rate for Adam optimizer.")
    p.add_argument("--dgi-scale", type=float, default=0.1,
                   help="Weight for DGI loss in total loss.")
    p.add_argument("--runs", type=int, default=1,
                   help="Number of independent runs to aggregate.")
    return p.parse_args()

# -----------------------------------------------------------------------------
# 8) Entry point
# -----------------------------------------------------------------------------

def main():
    start_time = time.time()   # 记录开始时间
    
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    for run_idx in range(1, args.runs + 1):
        print(f"\n===== Run {run_idx}/{args.runs} =====")
        run_outdir = args.outdir / f"run_{run_idx:02d}"
        train_single_run(
            edge_pairs_true_path=args.edge_pairs_true,
            edge_feats_true_path=args.edge_feats_true,
            edge_pairs_path=args.edge_pairs,
            edge_feats_path=args.edge_feats,
            node_emb_path=args.node_emb,
            pos_edges_path=args.pos_edges,
            neg_edges_path=args.neg_edges,
            outdir=run_outdir,
            epochs=args.epochs,
            hidden_dim=args.hidden_dim,
            heads=args.heads,
            sample_size=args.sample_size,
            lr=args.lr,
            dgi_scale=args.dgi_scale,
        )

    if args.runs > 1:
        aggregate_runs(args.outdir, args.runs, args.epochs)
        
    # 程序结束时打印消耗
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"\n✅ Total training time: {elapsed/60:.2f} minutes ({elapsed:.1f} seconds)")

if __name__ == "__main__":
    main()
