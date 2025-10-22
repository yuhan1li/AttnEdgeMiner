# AttnEdgeMiner

**AttnEdgeMiner** is a semi-supervised graph neural network (GNN) framework that leverages attention mechanisms to identify stable, high-confidence regulatory edges in biological networks, providing interpretable insights into gene interactions and signaling pathways.

---

## 📦 Requirements

- **R ≥ 4.2**  
- **Python ≥ 3.8**  
- **R packages**: `Seurat`, `dplyr`, `tidyverse`  
- **Python packages**: `torch`, `torch_geometric`, `node2vec`, `numpy`, `pandas`

---

## 🔹 Analysis Workflow

### Step 1: Run Network Construction (R)

```r
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

