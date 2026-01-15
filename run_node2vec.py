#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import time
import numpy as np
import torch
from torch_geometric.nn import Node2Vec


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train Node2Vec embeddings using PyTorch Geometric"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to edge list file. Format: node1<TAB>node2 (undirected)."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to save learned embeddings (TSV)."
    )
    parser.add_argument(
        "--embedding_dim",
        type=int,
        default=512,
        help="Embedding dimension (default: 512)."
    )
    parser.add_argument(
        "--walk_length",
        type=int,
        default=20,
        help="Random walk length (default: 20)."
    )
    parser.add_argument(
        "--context_size",
        type=int,
        default=10,
        help="Context window size (default: 10)."
    )
    parser.add_argument(
        "--walks_per_node",
        type=int,
        default=10,
        help="Number of walks per node (default: 10)."
    )
    parser.add_argument(
        "--p",
        type=float,
        default=0.5,
        help="Return parameter p (default: 0.5)."
    )
    parser.add_argument(
        "--q",
        type=float,
        default=1.0,
        help="In-out parameter q (default: 1.0)."
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=0.01,
        help="Learning rate (default: 0.01)."
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=200,
        help="Training epochs (default: 200)."
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=128,
        help="Batch size for Node2Vec loader (default: 128)."
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        default=4,
        help="Number of DataLoader workers (default: 4)."
    )
    parser.add_argument(
        "--gpu",
        type=int,
        default=0,
        help="GPU id to use (default: 0)."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Device setup
    if torch.cuda.is_available():
        device = torch.device(f"cuda:{args.gpu}")
        torch.cuda.set_device(device)
        print("Using GPU:", torch.cuda.get_device_name(args.gpu))
    else:
        device = torch.device("cpu")
        print("Using CPU")

    # Load edge list
    print("Loading edge list from:", args.input)
    try:
        edge_list = np.loadtxt(args.input, delimiter="\t", dtype=np.int64)
        edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
        edge_index = edge_index.to(device)
    except Exception as e:
        raise RuntimeError(f"Failed to load edge list: {e}")

    # Initialize Node2Vec model
    model = Node2Vec(
        edge_index=edge_index,
        embedding_dim=args.embedding_dim,
        walk_length=args.walk_length,
        context_size=args.context_size,
        walks_per_node=args.walks_per_node,
        num_negative_samples=1,
        p=args.p,
        q=args.q,
        sparse=True
    ).to(device)

    # Data loader and optimizer
    loader = model.loader(
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=args.num_workers
    )
    optimizer = torch.optim.SparseAdam(model.parameters(), lr=args.lr)

    # Training loop
    def train_one_epoch():
        model.train()
        total_loss = 0.0

        for pos_rw, neg_rw in loader:
            optimizer.zero_grad()
            loss = model.loss(
                pos_rw.to(device),
                neg_rw.to(device)
            )
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        return total_loss / len(loader)

    print("Start training Node2Vec")
    start_time = time.time()

    for epoch in range(1, args.epochs + 1):
        loss = train_one_epoch()

        if device.type == "cuda":
            allocated = torch.cuda.memory_allocated(device) / 1024 ** 2
            reserved = torch.cuda.memory_reserved(device) / 1024 ** 2
            peak = torch.cuda.max_memory_allocated(device) / 1024 ** 2

            print(
                f"Epoch {epoch:03d} | "
                f"Loss {loss:.4f} | "
                f"GPU alloc {allocated:.1f} MB | "
                f"reserved {reserved:.1f} MB | "
                f"peak {peak:.1f} MB"
            )
        else:
            print(f"Epoch {epoch:03d} | Loss {loss:.4f}")

    elapsed = time.time() - start_time
    print(f"Training finished in {elapsed:.2f} sec ({elapsed / 60:.2f} min)")

    if device.type == "cuda":
        peak_mem = torch.cuda.max_memory_allocated(device) / 1024 ** 2
        print(f"Peak GPU memory usage: {peak_mem:.1f} MB")

    # Save embeddings
    model.eval()
    embeddings = model.embedding.weight.detach().cpu().numpy()
    np.savetxt(args.output, embeddings, delimiter="\t", fmt="%.6f")
    print("Embeddings saved to:", args.output)


if __name__ == "__main__":
    main()

