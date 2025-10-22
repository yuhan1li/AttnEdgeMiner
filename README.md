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

The run_network_pipeline function constructs training data for graph-based models from single-cell RNA-seq data. It builds cell type-specific networks by integrating multiple data sources including prior knowledge (gene regulatory interactions, protein-protein interactions, and signaling pathways), highly expressed marker genes, and housekeeping genes. The function generates both positive edges (biologically validated gene pairs within each cell type) and negative edges (gene pairs with zero expression or random pairs lacking prior biological evidence). For each edge, it provides comprehensive information including node indices, edge features, and labels, creating a complete dataset optimized for downstream network learning or graph neural network training. The final output files are saved in the specified output_dir directory.

```r
# Load the run_network_pipeline function
# Source the required functions from the current directory
source("network_pipeline_function.R")

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
## 🔹 System requirements：

This model is developed and tested on a Linux GPU server with the following configurations: NVIDIA RTX A6000 GPUs and CUDA 11.7. The system uses NVIDIA driver version 550.107.02 (compatible with CUDA 12.4).
To facilitate environment setup and reproducibility, a complete environment.yml file is provided. Users can directly create the Python environment using this file, ensuring that all required dependencies are installed consistently for model execution and analysis.

## 🔹 创建python 环境
```
conda env create -f environment.yml
conda activate scDNS_python
```


### 2️⃣ Node Embedding with Node2Vec (Python)

run_node2vec.py generates low-dimensional embeddings for genes in the prior knowledge network. Each gene (node) is mapped to a 512-dimensional vector that captures both its local neighborhood structure and global network topology. These embeddings can then be used as input features for downstream tasks such as graph neural network training or edge prediction.

```
# Navigate to the output directory from the previous step
cd ./output_data/data_export/
# Place all Python script files in this folder for the next step execution

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
Use the --help argument to see all available input parameters and their descriptions.

### 3️⃣ Training Graph Attention + DGI (Python)

The train_gatv_2_dgi_multi.py script implements a sophisticated graph learning approach that combines GATv2 graph attention networks with Deep Graph Infomax (DGI). This integrated framework takes node embeddings and edge features generated in previous steps to perform edge classification, predicting potential regulatory relationships between gene pairs. The graph attention mechanism dynamically aggregates information from each node and its neighbors, while the self-supervised DGI module ensures the learned embeddings effectively capture the overall network structure. Supervised edge classification is optimized using positive and negative edge labels, with Focal Loss employed to handle class imbalance. The training outputs, including attention scores for each edge and refined node embeddings, are systematically saved in ./experiments/multi_run/ for subsequent network analysis and visualization.


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
Use the --help argument to see all available input parameters and their descriptions.
### 4️⃣ Extract and Merge Attention Scores (Python)

The combine_attention_scores.py script aggregates attention scores from multiple runs, integrates cell-type specific information, and computes comprehensive per-edge ranking metrics. It generates two key output files: a detailed per-edge score file (edge_scores_detailed.csv) containing raw attention scores from all runs, and a confidence summary file (edge_scores_confidence.csv) that provides consolidated rankings and confidence assessments for downstream biological interpretation.

```
python combine_attention_scores.py \
  --epoch 1000 \
  --root-dir ./experiments/multi_run \
  --celltype-file ./symbol_pair.txt \
  --output ./edge_scores_detailed.csv \
  --output-wide ./edge_scores_confidence.csv
```

### Final Attention Score Files

After completing the GATv2+DGI training and merging multiple runs, the final output files provide a comprehensive view of predicted gene-gene regulatory interactions.

- **Detailed per-edge scores (`edge_scores_confidence.csv`)**  
  This file contains per-edge attention scores from all training runs, along with computed ranking metrics within each cell type. Each row represents a gene pair (edge), and columns include:
  - `edge`: the gene pair, e.g., `AAAS_AGO2`
  - `KO_0` … `WT_6`: attention scores from different experimental conditions or replicates
  - `rank`: minimum rank across conditions
  - `final_rank`: overall rank across all cell types
  - `confidence`: high- or low-confidence classification based on cross-cell-type consistency

Example rows:

| edge       | KO_0   | KO_1   | ... | WT_6   | rank   | final_rank | confidence      |
|------------|--------|--------|-----|--------|--------|------------|----------------|
| AAAS_AGO2  | 123712 | 122726 | ... | 133692 | 117660 | 137408     | High_confidence |
| A1BG_AKT1  | 188970 | 189072 | ... | 198424 | 187237 | 195399     | Low_confidence  |





