# 🌟 AttnEdgeMiner

**AttnEdgeMiner** is a **semi-supervised graph neural network (GNN) framework** that leverages **attention mechanisms** to identify **stable, high-confidence regulatory edges** in biological networks. It provides interpretable insights into gene interactions and signaling pathways.  

---

## 📦 Requirements

### R Environment

- **R version ≥ 4.2**  
  - Recommended: **R 4.2.3** or later  
- **Required R packages**:  
  - `Seurat (v4)` – single-cell analysis  
  - `dplyr` – data manipulation  
  - `tidyverse` – collection of data science packages  
  - `Matrix` – sparse matrix operations  
  - `purrr` – functional programming  
  - `COSG` – marker gene identification  
  - `GSEABase` – handling gene sets (for housekeeping genes)

---

## 🔹 Analysis Workflow

The workflow consists of four main steps: **network construction**, **node embedding**, **GNN training**, and **attention score extraction**.

---

## 🔹 Example Data for Demonstration

To run the full workflow, we provide **two example datasets** for download. These datasets allow you to construct the Seurat object and the network object for analysis:

1. **Seurat object (`seurat_obj = iAT2_data`)**  
   - This is a Seurat V4 object containing processed single-cell data.  
   - Download link: [iAT2_data.rds](#)  

2. **scDNS output object (`iAT2_scDNSob`)**  
   - This object contains the network computed using scDNS for the same dataset, including MI and DREMI values.  
   - Download link: [iAT2_scDNSob.rds](#)  

> **Note:** Make sure to download both files before starting the analysis pipeline.

### 1️⃣ Network Construction (R)

run_network_pipeline constructs training data for graph-based models from single-cell RNA-seq data. It builds cell type–specific networks using prior knowledge (gene regulatory interactions, protein-protein interactions, and signaling pathways) and integrates highly expressed marker genes and housekeeping genes. The function generates positive edges (biologically relevant gene pairs within each cell type) and negative edges (gene pairs with zero expression or random pairs lacking prior evidence). For each edge, it provides node indices, edge features, and labels, producing a complete dataset ready for downstream network learning or GNN training. 最终生成的用于后续训练的文件保存在 output_dir 路径下面

```r
加载run_network_pipeline 函数
在当前路径下加载相关函数，需要
source network_pipeline_function.R

result_all <- run_network_pipeline(
  seurat_obj           = iAT2_data,
  group_col            = "WT_DSPKD",
  group_values         = c("WT", "KO"),
  network_obj          = iAT2_scDNSob@Network,
  random_network_obj   = iAT2_scDNSob@NEAModel[["DegreeData"]][["Network"]],
  housekeeping_gmt     = "./data/HSIAO_HOUSEKEEPING_GENES.v2025.1.Hs.gmt",
  output_dir           = "./output_data",
  run_clustering       = TRUE,
  dims                 = 1:20,
  resolution           = 0.5,
  top_n_marker         = 300,
  negative_sample_size = 50000
)
```

## 🔹 创建python 环境
```
conda env create -f environment.yml
conda activate scDNS_python
```


### 2️⃣ Node Embedding with Node2Vec (Python)

run_node2vec.py generates low-dimensional embeddings for genes in the prior knowledge network. Each gene (node) is mapped to a 512-dimensional vector that captures both its local neighborhood structure and global network topology. These embeddings can then be used as input features for downstream tasks such as graph neural network training or edge prediction.

```
进入第一步的输出文件路径下面进行进一步的基因embedding的获取
cd ./output_data/data_export/  并将所有的.py函数文件放在该文件夹下面， 方面下一步骤的执行

python run_node2vec.py \
  --input ./gene_map.txt \
  --output ./gene_embedding.txt \
  --embedding_dim 512 \
  --walk_length 20 \
  --context_size 10 \
  --walks_per_node 10 \
  --p 0.5 \
  --q 1.0 \
  --lr 0.01 \
  --epochs 200
```
Please use the argument --help to see all available input parameters.

### 3️⃣ Training Graph Attention + DGI (Python)

We use the train_gatv_2_dgi_multi.py script for further analysis. In this step, the gene regulatory network is trained using a GATv2 graph attention network combined with Deep Graph Infomax (DGI). The script takes the node embeddings and edge feature files generated in the previous steps to perform edge classification, predicting potential regulatory relationships between gene pairs. The graph attention mechanism aggregates information from each node and its neighbors, while the self-supervised DGI module ensures the embeddings capture the overall network structure. Supervised edge classification is optimized with positive and negative edge labels, using Focal Loss to handle class imbalance. The outputs, including attention scores for each edge and updated node embeddings, are saved in ./experiments/multi_run/ for downstream network analysis and visualization.

```
python train_gatv_2_dgi_multi.py \
  --edge-pairs-true ./index_pairs_true.txt \
  --edge-feats-true ./edge_feature_true.txt \
  --edge-pairs ./index_pairs_all.txt \
  --edge-feats ./edge_feature_all.txt \
  --node-emb ./gene_embedding.txt \
  --pos-edges ./Positive_edge.txt \
  --neg-edges ./Negative_edge.txt \
  --outdir ./experiments/multi_run \
  --epochs 1000 \
  --runs 10 \
  --dgi-scale 0.1
```
Please use the argument --help to see all available input parameters.
### 4️⃣ Extract and Merge Attention Scores (Python)
```
python combine_attention_scores.py \
  --epoch 1000 \
  --root-dir ./experiments/multi_run \
  --celltype-file ./symbol_pair.txt \
  --output ./attention_layer1_epoch1000_combined_rank.csv \
  --output-wide ./attention_layer1_epoch1000_wide_confidence.csv
```





