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

```r
加载run_network_pipeline 函数
source network_pipeline_function.R

result_all <- run_network_pipeline(
  seurat_obj           = iAT2_data,
  group_col            = "WT_DSPKD",
  group_values         = c("WT", "KO"),
  network_obj          = iAT2_scDNSob@Network,
  random_network_obj   = iAT2_scDNSob@NEAModel[["DegreeData"]][["Network"]],
  housekeeping_gmt     = "../gene_embedding/HSIAO_HOUSEKEEPING_GENES.v2025.1.Hs.gmt",
  output_dir           = "./data",
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
```
cd ./data/data_export/

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

### 3️⃣ Training Graph Attention + DGI (Python)
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

### 4️⃣ Extract and Merge Attention Scores (Python)
```
python combine_attention_scores.py \
  --epoch 1000 \
  --root-dir ./experiments/multi_run \
  --celltype-file ./symbol_pair.txt \
  --output ./attention_layer1_epoch1000_combined_rank.csv \
  --output-wide ./attention_layer1_epoch1000_wide_confidence.csv
```





