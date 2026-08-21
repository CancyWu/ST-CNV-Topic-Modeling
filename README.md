# Interpretable Spatial Cancer Clone Analysis Using Slide-DNA-seq

## Initial Methods and Results Report

## 1. Objective

To identify spatial cancer clone-associated CNA patterns and connect the resulting spatial groups to interpretable signed genomic features and regions.

## 2. Key outcomes

### 2.1 Final model output

The final model learned four overlapping CNV topics from 30,392 Slide-DNA-seq beads and produced three spatial CNA groups by downstream clustering of the topic proportions. The model was trained without clone labels.

The three-group partition was reproducible across independent training runs (mean pairwise ARI = 0.933). The groups were spatially organized and showed distinct mean CNA profiles.

![Final spatial groups and dominant CNV topics](../Step12_sparsity_optimization/output_candidate/selected_l1tv010/K4_spatial_clones_and_topics.png)

### 2.2 Biological interpretation

The signed loading matrix linked each topic to gain-like and loss-like genomic regions. The strongest topic-associated patterns included large regions on chromosomes 8, 15, and 20. Adjacent informative 1 Mb bins were aggregated into continuous CNA regions for interpretation.

<img width="3181" height="2347" alt="topic_region_loading_shares" src="https://github.com/user-attachments/assets/5f382957-3aff-4ccc-8834-fc0698810e13" />

### 2.3 External validation

Matched scWGS independently supported the broad biological identities of the three spatial groups:

- Spatial group 1 matched a tumour CNA reference;
- Spatial group 2 matched a different tumour CNA reference;
- Spatial group 3 matched the normal CNA reference.

The corresponding profile correlations were 0.767, 0.635, and 0.862. Bead-level agreement was moderate across the complete section but substantially stronger in high-confidence spatial-interior regions, indicating that most uncertainty was concentrated near spatial boundaries and transition areas.

<img width="2640" height="1200" alt="spatial_model_vs_scwgs_projection" src="https://github.com/user-attachments/assets/862378ae-2d07-44d9-ad91-8c6868d2a0ed" />

### 2.4 Summary of initial conclusions

The initial analysis supports three main conclusions:

1. The CNV features contain sufficient information to recover the main spatially organized groups.
2. The true spatial graph improves reproducibility and feature stability rather than simply creating visually smooth regions.
3. The resulting groups have interpretable signed CNA profiles and receive independent support from matched scWGS.

## 3. Model workflow

### 3.1 End-to-end training

```text
CNV features X + spatial graph
              ↓
         GAT encoder
              ↓
Topic proportions Z = softmax(GAT output)
              ↓
CNV reconstruction: X_hat = baseline + ZB
              ↑
Signed topic-feature loading B
              ↓
Reconstruction loss end-to-end updates GAT and B
```

### 3.2 Spatial-group assignment

```text
Final topic proportions Z
              ↓
       K-means clustering
              ↓
Three spatial CNA groups
              ↓
Group-level CNA profiles and biological interpretation
```

The topic model was trained without clone labels. K-means was applied only after end-to-end training as a downstream discretization of the continuous topic representation.

### 3.3 Validation overview

Validation focused on four aspects:

- **Model fit and interpretability:** reconstruction, loading sparsity, and region-level CNA interpretation;
- **Robustness:** stability of spatial groups and signed genomic features across training seeds;
- **Contribution of spatial information:** comparison with CNV-only and shuffled-graph controls;
- **External validation:** CNA-profile matching and bead-level projection using matched scWGS, including confidence and spatial-boundary analyses.

## 4. Methods

### 4.1 Input and preprocessing

The raw input consisted of bead-by-1 Mb genomic-bin DNA fragment counts and bead spatial coordinates. Beads outside the analysed puck radius or with no more than 100 nuclear fragments were removed. Counts were aggregated over the 50 nearest spatial beads and normalized by the mean genome-wide signal of each spatial unit.

The final input contained:

- 30,392 spatial units;
- 2,704 autosomal 1 Mb genomic bins;
- normalized relative coverage for every bead-bin pair;
- two-dimensional spatial coordinates.

A signal near 1 indicates coverage close to the bead-specific genome-wide baseline. Values above or below 1 indicate relative gain-like or loss-like signals. These values are relative coverage measurements rather than absolute integer copy numbers.

### 4.2 GAT spatial encoder

A spatial 10-nearest-neighbour graph was constructed from bead coordinates. A two-layer GATv2 encoder propagated information between spatial neighbours and produced a spatially informed representation for every bead. The encoder used a 128-dimensional hidden layer and a 16-dimensional latent representation.

Only the spatial graph was used in the primary model; no CN-profile similarity graph was added. Genomic features were standardized and clipped before entering the encoder to reduce the influence of extreme bins.

### 4.3 Topic proportions and signed decoder

The final model used four CNV topics. Topic proportions were calculated as:

```text
Z = softmax(X_GAT)
```

Each row of `Z` sums to one. Multiple topics can therefore contribute to the same bead.

The signed topic-feature loading matrix `B` has dimensions `4 topics × 2,704 genomic bins` and links topics to genomic features. Positive and negative values represent gain-like and loss-like contributions. The normalized CNV input was reconstructed as:

```text
X_hat = baseline + ZB
```

Reconstruction loss was backpropagated through `B`, `Z`, and the GAT encoder, allowing spatial representation and genomic interpretation to be learned jointly. Mean-squared error was used because the input consists of continuous normalized coverage rather than counts.

The loss additionally included L1 regularization for loading concentration, genomic total variation for locally compatible loadings, and small topic-entropy and topic-balance terms. The selected model used `L1 = 0.10` and `TV = 0.01`.

### 4.4 Region-level interpretation

Interpretation was performed at two levels:

- **Bin level:** signed loadings identify influential 1 Mb genomic bins.
- **Region level:** adjacent bins with concordant signed loadings are merged into continuous CNA regions.

Region aggregation was post-hoc and did not modify the trained model. The primary analysis used the topic-specific 85th percentile of absolute loading, required at least two selected bins, and allowed one weak intervening bin. The 80th- and 90th-percentile results were retained as sensitivity analyses.

### 4.5 scWGS validation design

The published Slide-DNA-seq spatial labels were not used as ground truth. Instead, matched scWGS provided an independent CNA reference library containing one normal and five tumour profiles. Each spatial group was matched to the scWGS library by genomic-profile correlation. The best-matching normal reference and two tumour references were then projected onto the Slide-DNA-seq beads.

## 5. Detailed results and validation

### 5.1 Reconstruction, sparsity, and training stability

The final model achieved a reconstruction RMSE of 0.3689, compared with 0.3793 for the feature-mean baseline. Sparsity optimization increased mean loading Gini from 0.474 to 0.669 and reduced the effective number of contributing bins from approximately 1,805 to 1,144. The top 10% of bins accounted for 47.4% of total absolute loading, compared with 33.3% before optimization.

Across three training seeds, the three-group partitions had a mean pairwise ARI of 0.933 and a minimum ARI of 0.931.

### 5.2 Stability of genomic interpretation

After matching topics across seeds:

- mean selected-bin Jaccard overlap was 0.905;
- direction consistency was 1.000;
- 225 stable consensus regions were detected.

Thus, both the locations and gain/loss directions of influential genomic features were reproducible across training runs.

### 5.3 Contribution of spatial information

| Metric | Spatial | CNV-only | Shuffled |
|---|---:|---:|---:|
| Across-seed cluster ARI | **0.933** | 0.608 | 0.924 |
| Top-loading feature overlap | **0.905** | 0.845 | 0.606 |
| Raw-CNV PCA silhouette | **0.171** | 0.165 | 0.011 |
| Clone-profile distance | **0.1675** | 0.1669 | 0.0864 |
| Reconstruction RMSE | **0.3689** | 0.3692 | 0.3779 |

Spatial and CNV-only partitions remained similar (ARI = 0.867), indicating that the main group identities were present in the CNV features. The true spatial graph improved cross-seed partition and feature stability. The shuffled graph produced stable but poorly separated groups, showing that arbitrary graph smoothing could not replace genuine spatial adjacency.

<img width="2860" height="1760" alt="ablation_metric_dashboard" src="https://github.com/user-attachments/assets/b29df4b8-5b95-4829-9e25-257096e0b785" />

### 5.4 scWGS profile matching

The two scWGS batches contained 2,341 raw cells and 2,274 QC-passing cells. After excluding 357 doublet candidates, 1,917 cells formed the normal-plus-tumour reference library. The two assays were aligned over 2,435 shared high-quality 1 Mb bins.

| Spatial group | Best-matching scWGS reference | Pearson correlation |
|---|---|---:|
| Spatial group 1 | Tumour 5 | 0.767 |
| Spatial group 2 | Tumour 1 | 0.635 |
| Spatial group 3 | Normal | 0.862 |

<img width="1760" height="836" alt="clone_profile_correlation_heatmap" src="https://github.com/user-attachments/assets/31f628c3-e0d5-4cfb-8cd6-17ebb09cc5b1" />

The three-reference bead projection produced ARI = 0.589, NMI = 0.597, and AMI = 0.597 across the complete section.

### 5.5 Projection confidence and spatial boundaries

Whole-section disagreement was concentrated in uncertain and boundary beads.

| Bead subset | ARI | NMI | AMI | Direct agreement |
|---|---:|---:|---:|---:|
| All beads | 0.589 | 0.597 | 0.597 | 84.9% |
| Spatial-interior beads | 0.763 | 0.755 | 0.755 | 91.9% |
| Spatial-boundary beads | 0.298 | 0.326 | 0.326 | 65.0% |
| Interior and high-confidence | 0.893 | 0.876 | 0.876 | 96.8% |

Spatial group 3 was projected as normal for 98.7% of its beads. Spatial groups 1 and 2 matched their expected tumour references for 78.6% and 83.0% of beads. These results support the core group identities while identifying spatial transitions as the main source of uncertainty.

<img width="2640" height="1032" alt="projection_confidence_and_confusion" src="https://github.com/user-attachments/assets/e0b68fbd-7ee8-4d0f-995c-944120f57dfc" />

### 5.6 Region-level topic interpretation

At the primary threshold, aggregated regions accounted for 45.5%, 52.8%, 51.2%, and 57.8% of total absolute loading in Topics 1–4.

- Topic 1: loss-like chr8:58–113 Mb and chr20 regions;
- Topic 2: loss-like chr8:47–146 Mb and weaker gain-like chr15 regions;
- Topic 3: gain-like chr8:67–111 Mb and chr20 regions;
- Topic 4: gain-like chr8:47–146 Mb and loss-like chr15:48–82 Mb.

The opposite signed chr8 patterns show that the model captures gain-like and loss-like genomic programs rather than only non-negative feature importance.

## 6. Initial conclusion

The current model provides an end-to-end spatial topic representation of normalized Slide-DNA-seq coverage. It identifies three stable spatial groups, connects topics to signed genomic regions, and receives independent support from matched scWGS. The main CNA identities are present in the CNV data, while genuine spatial information improves reproducibility and feature stability. External agreement is strongest in high-confidence interior regions, with uncertainty concentrated at spatial boundaries.

These results support the feasibility of interpretable spatial clone analysis, although neither the inferred spatial groups nor the scWGS references should be treated as absolute ground truth.
