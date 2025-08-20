# Statistical Inference for Subgraph Frequencies of Exchangeable Hyperedge Models

## Overview

In statistical network analysis, models for binary adjacency matrices satisfying vertex exchangeability are commonly used. However, such models may fail to capture key features of the data-generating process when interactions, rather than nodes, are fundamental units. We study statistical inference for subgraph counts under an exchangeable hyperedge model.  We introduce several classes of subgraph statistics for hypergraphs and develop inferential tools for subgraph frequencies that account for edge multiplicity. We show that a subclass of these subgraph statistics is robust to the deletion of low-degree nodes, enabling inference in settings where low-degree nodes are more likely to be missing. We also examine a more traditional notion of subgraph frequency that ignores multiplicity, showing that while inference based on limiting distributions is feasible in some cases, a non-degenerate limiting distribution may not exist in others. Empirically, we assess our methods through simulations and newly collected real-world hypergraph data on academic and movie collaborations, where our inferential tools outperform traditional approaches based on binary adjacency matrices.

## Main Article

The full article is available a https://arxiv.org/abs/2508.13258.

## Repository Structure

```
hypergraph_subgraph/
│
├── figures/ # All figures for the paper
│
├── real_data/
│ ├── data_collection/ # Real-world hypergraph datasets (.RData)
│ ├── hypergraph_testing/ # Tests based on hypergraph framework
│ └── sparse_graphon_testing/# Tests based on sparse graphon models
│
├── simulation/ # Simulation experiments
│
└── README.md # This file
```


## Requirements

The code is written in **R** and requires the following packages:

```r
install.packages(c(
  "igraph",
  "rje",
  "ggplot2",
  "parallel",
  "RSpectra"
))
```

## Reproducing the Results

1. Clone the repository:

```bash
git clone https://github.com/ayoushmanb/hypergraph_subgraph.git
cd hypergraph_subgraph
```

## Run simulations:

```r
source("simulation/sim-colorless-triangles-run.R")
source("simulation/sim-type2-clustering-coeff-run.R")
source("simulation/sim-type2-triangle-run.R")
source("simulation/sim-type2-twostar-run.R")
source("simulation/sim-type3-clustering-coeff-run.R")
source("simulation/sim-type3-triangle-run.R")
```

## Run real data analysis:

```r
source("real_data/hypergraph_testing/real-data-twostar2-test.R")
source("real_data/hypergraph_testing/real-data-type2-clustering-test.R")
source("real_data/sparse_graphon_testing/real-data-sparse-graphon-clus-coeff.R")
source("real_data/sparse_graphon_testing/real-data-sparse-graphon-two-stars.R")
```

## Data

The real_data/data_collection/ folder contains preprocessed .RData files for:

  - EMI, JASA, JCGS, and Movie collaboration datasets
  - Each dataset has *_train.RData, *_test.RData, and *_all.RData versions.

## Citation
If you use this code or data, please cite:

```mathematica
@article{bhattacharya2025statistical,
  title={{Statistical Inference for Subgraph Frequencies of Exchangeable Hyperedge Models}},
  author={Bhattacharya, Ayoushman and Chakraborty, Nilanjan and Lunde, Robert},
  journal={arXiv preprint},
  volume={arXiv:2508.13258},
  year={2025},
  month={Aug},
  note={\url{https://doi.org/10.48550/arXiv.2508.13258}}
}
```
