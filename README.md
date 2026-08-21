# Interpretable Spatial Cancer Clone Analysis Using Slide-DNA-seq

## Initial Methods and Results Report

## 1. Objective

To develop an interpretable spatial topic model that identifies cancer clone-associated CNA patterns from Slide-DNA-seq while linking each topic to signed genomic features and regions.

## 2. Methods

### 2.1 Model and analysis overview

#### End-to-end training framework

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

#### Clone-group assignment after training

```text
Final topic proportions Z
              ↓
       K-means clustering
              ↓
Three spatial CNA groups
              ↓
Group-level CNA profiles and biological interpretation
```

The topic model is trained without clone labels. K-means is applied only after the end-to-end topic model has been trained and is therefore a downstream discretization of the continuous topic representation.

#### Clarification of each component

1. **CNV features \(X\).** Each spatial unit is represented by normalized relative DNA coverage across 2,704 autosomal 1 Mb genomic bins. The values are continuous relative signals, not absolute integer copy numbers.

2. **Spatial graph.** Beads are connected through a spatial 10-nearest-neighbour graph constructed from their physical coordinates. No CN-profile similarity graph is included in the primary model.

3. **GAT encoder.** The GAT propagates information between spatial neighbours and produces a spatially informed representation for every bead.

4. **Topic proportions \(Z\).** A softmax transformation produces four non-negative topic proportions per bead. Each row of \(Z\) sums to one, allowing topics to overlap within a bead rather than assigning every bead to a single topic.

5. **Signed loading matrix \(B\).** Each row of \(B\) describes the contribution of genomic bins to one topic. Positive and negative values represent gain-like and loss-like contributions, respectively.

6. **CNV reconstruction.** The model reconstructs the original normalized CNV input as a feature baseline plus the product of bead-topic proportions and topic-feature loadings.

7. **End-to-end optimization.** Reconstruction loss is backpropagated through $$\widehat{X}$$, \(B\), \(Z\), and the GAT encoder. Thus, the spatial representation and genomic interpretation are learned jointly.

8. **Clone-group assignment.** K-means with three groups is applied to final \(Z\) for a discrete spatial partition. The continuous topic proportions remain available for analysing mixed or overlapping topic states.

#### Validation analyses completed

Validation focused on four aspects:

- **Model fit and interpretability:** reconstruction, loading sparsity, and region-level CNA interpretation;
- **Robustness:** stability of spatial groups and signed genomic features across training seeds;
- **Contribution of spatial information:** comparison with CNV-only and shuffled-graph controls;
- **External validation:** CNA-profile matching and bead-level projection using matched scWGS, including confidence and spatial-boundary analyses.

### 2.2 Slide-DNA-seq input and preprocessing

The raw input consisted of bead-by-1 Mb genomic-bin DNA fragment counts and bead spatial coordinates. Beads outside the analysed puck radius or with no more than 100 nuclear fragments were removed. Counts were aggregated over the 50 nearest spatial beads and normalized by the mean genome-wide signal of each spatial unit.

After quality control, the model input contained:

- 30,392 spatial units (beads);
- 2,704 autosomal 1 Mb genomic bins;
- a normalized relative coverage value for every bead-bin pair;
- two-dimensional spatial coordinates for every bead.

A value near 1 represents coverage close to the bead-specific genome-wide baseline. Values above or below 1 indicate relative gain-like or loss-like signals. These features are relative coverage measurements rather than absolute integer copy numbers.

### 2.3 Spatial encoding

A spatial 10-nearest-neighbour graph was constructed from bead coordinates. A two-layer GATv2 encoder was used to propagate information between spatial neighbours:

```text
Normalized CNV features + spatial 10-NN graph
                         ↓
                 GATv2 spatial encoder
                         ↓
                   spatial representation
```

The encoder used a 128-dimensional hidden layer followed by a 16-dimensional latent representation. Only the spatial graph was used; no CN-profile similarity graph was added. Before encoding, each genomic feature was standardized and clipped to reduce the influence of extreme bins.

### 2.4 Topic representation and signed feature decoder

The final model used four CNV topics. Topic proportions for each bead were obtained by applying a softmax function to the GAT output:

$$
Z = \text{softmax}(X_{\text{GAT}})
$$

where $Z$ has dimensions spatial unit × topic and each row sums to one.

To make the topics interpretable, a signed topic-feature loading matrix $B$ was learned:

$$
B \in \mathbb{R}^{K \times P}
$$

where $K=4$ topics and $P=2{,}704$ genomic bins. Positive and negative loadings represent gain-like and loss-like contributions, respectively. The normalized CNV input was reconstructed as

$$
\hat{X} = \mu + ZB
$$

where $\mu$ is the fixed feature-wise baseline.

The GAT encoder, topic proportions, and signed loading matrix were optimized end-to-end using a single reconstruction objective. Because the input is continuous normalized coverage rather than count data, mean-squared reconstruction error was used instead of a Gamma-Poisson likelihood.

The training objective included:

- CNV reconstruction loss;
- L1 regularization on \(B\) to concentrate loading mass;
- genomic total-variation regularization to encourage neighbouring bins to have compatible loadings;
- small topic-entropy and topic-balance terms.

The selected model used `L1 = 0.10` and `TV = 0.01`. Three downstream spatial groups were derived by applying K-means clustering to the final topic proportions. Clone labels were not used during model training.

### 2.5 Interpretation of genomic features

Interpretation was performed at two levels:

1. **Bin level:** positive and negative topic loadings identify influential 1 Mb genomic bins.
2. **Region level:** adjacent bins with concordant signed loadings were merged into continuous CNA regions.

Region aggregation was a post-hoc interpretation step and did not alter the trained model. The primary analysis selected bins above the topic-specific 85th percentile of absolute loading, required at least two selected bins per region, and allowed one weak intervening bin. Results at the 80th and 90th percentiles were retained as sensitivity analyses.

### 2.6 External scWGS validation design

The model was evaluated without treating the published Slide-DNA-seq spatial partition as ground truth.

For scWGS validation, one normal and five tumour CNA reference profiles were constructed from matched scWGS cells. The three spatial-model groups were matched to this reference library by genomic-profile correlation. The best-matching normal reference and two tumour references were then projected onto the Slide-DNA-seq beads.

## 3. Results

### 3.1 Reconstruction and sparsity

The final four-topic model achieved a reconstruction RMSE of 0.3689, compared with 0.3793 for the feature-mean baseline. Increasing L1 regularization from the initial setting improved loading concentration with little reconstruction penalty:

- mean loading Gini increased from 0.474 to 0.669;
- the effective number of contributing bins decreased from approximately 1,805 to 1,144;
- the top 10% of bins accounted for 47.4% of total absolute loading, compared with 33.3% before sparsity optimization.

Across three joint-training seeds, the three-group partitions had a mean pairwise ARI of 0.933 and a minimum ARI of 0.931. Thus, the selected regularization improved interpretability while preserving reconstruction and partition stability.

### 3.2 Stability of feature interpretation

Topics were aligned across training seeds using Hungarian matching based on topic-proportion and loading correlations. Important bins were defined as the top 15% by absolute loading within each topic and seed.

- Mean selected-bin Jaccard overlap across seeds: 0.905.
- Direction consistency of matched selected bins: 1.000.
- Stable consensus regions detected across seeds: 225.

These results indicate that both the identity and gain/loss direction of influential genomic features were reproducible across independent training runs.

### 3.3 Contribution of spatial information

The spatial ablation compared the true spatial graph with CNV-only/self-only encoding and a shuffled graph.

| Metric | Spatial | CNV-only | Shuffled |
|---|---:|---:|---:|
| Across-seed cluster ARI | **0.933** | 0.608 | 0.924 |
| Top-loading feature overlap | **0.905** | 0.845 | 0.606 |
| Raw-CNV PCA silhouette | **0.171** | 0.165 | 0.011 |
| Clone-profile distance | **0.1675** | 0.1669 | 0.0864 |
| Reconstruction RMSE | **0.3689** | 0.3692 | 0.3779 |

The Spatial and CNV-only representative partitions remained similar (ARI = 0.867), showing that the main clone signal was already present in the CNV features. However, the true spatial graph substantially improved across-seed partition stability and feature-loading stability. The shuffled graph produced stable but poorly separated groups, demonstrating that arbitrary graph smoothing could not replace genuine spatial adjacency.

<img width="2860" height="1760" alt="ablation_metric_dashboard" src="https://github.com/user-attachments/assets/b29df4b8-5b95-4829-9e25-257096e0b785" />

### 3.4 External validation using matched scWGS

The two scWGS batches contained 2,341 raw cells and 2,274 QC-passing cells. After excluding 357 doublet candidates, 1,917 cells formed a reference library containing one normal group and five tumour groups. The Slide-DNA-seq and scWGS profiles were aligned over 2,435 shared high-quality 1 Mb bins.

The three spatial groups matched the scWGS references as follows:

| Spatial group | Best-matching scWGS reference | Pearson correlation |
|---|---|---:|
| Spatial group 1 | Tumour 5 | 0.767 |
| Spatial group 2 | Tumour 1 | 0.635 |
| Spatial group 3 | Normal | 0.862 |

The result supports one normal-like spatial group and two tumour-like groups with distinct CNA profiles. The three matched references were projected onto individual beads, producing overall concordance of ARI = 0.589, NMI = 0.597, and AMI = 0.597.

<img width="1760" height="836" alt="clone_profile_correlation_heatmap" src="https://github.com/user-attachments/assets/31f628c3-e0d5-4cfb-8cd6-17ebb09cc5b1" />

<img width="2640" height="1200" alt="spatial_model_vs_scwgs_projection" src="https://github.com/user-attachments/assets/862378ae-2d07-44d9-ad91-8c6868d2a0ed" />

### 3.5 Confidence and spatial-boundary diagnostics

The moderate whole-section concordance was concentrated in uncertain and boundary beads. Projection confidence was defined as the correlation difference between the best and second-best scWGS references.

| Bead subset | Beads | ARI | NMI | AMI | Direct agreement |
|---|---:|---:|---:|---:|---:|
| All beads | 30,392 | 0.589 | 0.597 | 0.597 | 84.9% |
| Spatial-interior beads | 22,452 | 0.763 | 0.755 | 0.755 | 91.9% |
| Spatial-boundary beads | 7,940 | 0.298 | 0.326 | 0.326 | 65.0% |
| Interior and high-confidence | 11,903 | 0.893 | 0.876 | 0.876 | 96.8% |

Spatial group 3 was projected as normal for 98.7% of its beads. Spatial groups 1 and 2 matched their expected tumour references for 78.6% and 83.0% of their beads, respectively. These results suggest that the core regions and mean CNA identities are strongly supported, while uncertainty is concentrated near spatial transitions and in beads with similar correlations to multiple references.

<img width="2640" height="1032" alt="projection_confidence_and_confusion" src="https://github.com/user-attachments/assets/e0b68fbd-7ee8-4d0f-995c-944120f57dfc" />


### 3.6 Region-level topic interpretation

Aggregating adjacent signed bins produced concise region-level interpretations. At the primary 85th-percentile threshold, identified regions accounted for 45.5%, 52.8%, 51.2%, and 57.8% of the total absolute loading in Topics 1–4, respectively.

The dominant patterns included:

- Topic 1: loss-like chr8:58–113 Mb and chr20 regions;
- Topic 2: a dominant loss-like chr8:47–146 Mb region and weaker gain-like chr15 regions;
- Topic 3: gain-like chr8:67–111 Mb and chr20 regions;
- Topic 4: a dominant gain-like chr8:47–146 Mb region and loss-like chr15:48–82 Mb.

<img width="3181" height="2347" alt="topic_region_loading_shares" src="https://github.com/user-attachments/assets/5f382957-3aff-4ccc-8834-fc0698810e13" />

The opposite signed chr8 patterns across topics demonstrate that the loading matrix can distinguish gain-like and loss-like genomic programs rather than returning only non-negative feature weights.

## 4. Initial conclusion

The current model provides an end-to-end spatial topic representation of normalized Slide-DNA-seq coverage. It identifies three stable spatial groups, connects latent topics to signed genomic features, and produces region-level CNA interpretations. The main group identities are driven by the CNV data, while the true spatial graph improves reproducibility and feature stability. Matched scWGS independently supports one normal-like group and two distinct tumour-like groups, with the strongest agreement in high-confidence interior regions.

The present results support the feasibility of interpretable spatial clone analysis, but the inferred groups and scWGS references should not be treated as absolute ground truth. Further work may refine uncertainty handling, region definitions, and validation on additional datasets.
