#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import torch
import numpy as np
import time
from torch_geometric.nn import Node2Vec


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train Node2Vec embeddings using PyTorch Geometric"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to edge list file. Format: node1<TAB>node2 (each line represents an undirected edge)."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to save the learned embedding file (tab-separated values)."
    )
    parser.add_argument(
        "--embedding_dim", type=int, default=512,
        help="Size of the node embedding vector (default: 512)"
    )
    parser.add_argument(
        "--walk_length", type=int, default=20,
        help="Length of each random walk (default: 20)"
    )
    parser.add_argument(
        "--context_size", type=int, default=10,
        help="Context window size for the Skip-Gram model (default: 10)"
    )
    parser.add_argument(
        "--walks_per_node", type=int, default=10,
        help="How many random walks to start per node (default: 10)"
    )
    parser.add_argument(
        "--p", type=float, default=0.5,
        help="Node2Vec return parameter p (lower -> more likely to go back to previous node, default: 0.5)"
    )
    parser.add_argument(
        "--q", type=float, default=1.0,
        help="Node2Vec in-out parameter q (higher -> more exploratory search, default: 1.0)"
    )
    parser.add_argument(
        "--lr", type=float, default=0.01,
        help="Learning rate for SparseAdam optimizer (default: 0.01)"
    )
    parser.add_argument(
        "--epochs", type=int, default=200,
        help="Number of training epochs (default: 200)"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Check GPU
    if torch.cuda.is_available():
        device = torch.device('cuda')
        print("✅ Using GPU:", torch.cuda.get_device_name(0))
    else:
        device = torch.device('cpu')
        print("⚙️ Using CPU")

    # Load edge list
    print("📂 Loading edge list from:", args.input)
    try:
        edge_index = np.loadtxt(args.input, delimiter='\t', dtype=int).tolist()
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_index = edge_index.to(device)
    except Exception as e:
        raise ValueError(f"❌ Failed to load edge file: {e}")

    # Initialize Node2Vec model
    node2vec = Node2Vec(
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

    # Optimizer and dataloader
    loader = node2vec.loader(batch_size=128, shuffle=True, num_workers=4)
    optimizer = torch.optim.SparseAdam(list(node2vec.parameters()), lr=args.lr)

    # Training function
    def train():
        node2vec.train()
        total_loss = 0
        for pos_rw, neg_rw in loader:
            optimizer.zero_grad()
            loss = node2vec.loss(pos_rw.to(device), neg_rw.to(device))
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        return total_loss / len(loader)

    # Start training
    print("🚀 Training Node2Vec...")
    start = time.time()
    for epoch in range(1, args.epochs + 1):
        loss = train()
        print(f"Epoch {epoch:03d} | Loss: {loss:.4f}")

    elapsed = time.time() - start
    print(f"⏱️ Finished! Time cost: {elapsed:.2f} sec ({elapsed/60:.2f} min)")

    # Save embedding
    node2vec.eval()
    embeddings = node2vec.embedding.weight.detach().cpu().numpy()
    np.savetxt(args.output, embeddings, delimiter='\t', fmt="%.6f")
    print(f"✅ Embeddings saved to: {args.output}")


if __name__ == "__main__":
    main()
