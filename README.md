# FlyWire Qualification Challenge
**Largest Shared Neuronal Circuit Across Connectome Datasets**

Created by: Samira Salimiyan | PhD student in Biomedical Sciences | June 2026

---

## Overview

This repository contains the solution to the FlyWire Academy Qualification Challenge: identifying the largest neuronal circuit shared across at least three of the five FlyWire connectome datasets.

**Datasets selected:** FAFB (Female Adult Brain) + BANC (Female Brain + Nerve Cord) + MCNS (Male Full CNS)

**Approach:** Degree signature-based node correspondence followed by signature-level edge intersection to identify conserved directed connectivity motifs across independently reconstructed connectomes.

---

## Repository Contents

| File | Description |
|------|-------------|
| `flywire_circuit_submission.csv` | Solution file: 3 columns for the 3 selected datasets (neuron IDs), N rows (matched neurons) |
| `flywire_analysis_clean.R` | Full documented R analysis pipeline |
| `flywire_report.docx` | Scientific summary with circuit visualization and biological interpretation |
| `README.md` | This file |

---

## Methods Summary

### Dataset Selection
Five unweighted directed graphs were loaded from provided edge lists. Degree distributions were computed and visualized for all five datasets. FAFB, BANC, and MCNS were selected based on overlapping degree profiles (median total degree: 40, 24, and 46 respectively), making them the most compatible trio for connectivity-based matching. MANC (median 335) and MAOL (median 183) were excluded due to incompatible degree ranges reflecting their specialized nerve cord and optic lobe wiring.

### Node Correspondence via Degree Signatures
A degree signature was defined as the (in-degree, out-degree) pair for each neuron. Neurons were filtered to total degree 3–200 to capture the shared connectivity zone. Signatures present in all three datasets were identified using set intersection, yielding 6,419 shared signatures across FAFB, BANC, and MCNS.

### Circuit Identification
Edges from induced subgraphs of candidate neurons were translated to signature space, replacing neuron IDs with degree signatures to enable cross-dataset comparison. The intersection of signature-level edges across all three datasets defined conserved directed connectivity motifs (37,652 shared edges). A directed graph was constructed from these edges; the largest weakly connected component was extracted as the candidate circuit.

### Node Mapping
One representative neuron per degree-signature class per dataset was selected for the final correspondence table. This is valid because neurons sharing a degree signature are structurally interchangeable within the isomorphism framework. The final circuit contains N matched neurons across all three datasets.

---

## Key Results

| Metric | Value |
|--------|-------|
| Shared degree signatures | 6,419 |
| Shared signature-level edges | 37,652 |
| Circuit size (N neurons) | 2942 |
| Datasets | FAFB + BANC + MCNS |
| Connectivity | Weakly connected (largest component) |

---

## Limitations and Discussion

Degree signature matching is a necessary but not sufficient condition for strict graph isomorphism. Two neurons may share identical (in-degree, out-degree) pairs while connecting to different cell types. Full biological isomorphism requires cell-type annotation as the matching key, which is the gold standard used in published cross-connectome comparisons (Schlegel et al., 2023). The current solution identifies structurally conserved motifs at the degree level; verification against Codex cell type metadata would strengthen the correspondence.

---

## Requirements

```r
R >= 4.0
tidyverse
igraph
data.table
```

## Usage

```r
setwd("path/to/your/data/directory")
source("flywire_analysis_clean.R")
# Outputs saved as output_01 through output_11 in working directory
```

---

## References

- Dorkenwald S, McKellar CE, Macrina T, Kemnitz N, Lee K, Lu R, Wu J, Popovych S, Mitchell E, Nehoran B, Jia Z, Bae JA, Mu S, Ih D, Castro M, Ogedengbe O, Halageri A, Kuehner K, Sterling AR, Ashwood Z, Zung J, Brittain D, Collman F, Schneider-Mizell C, Jordan C, Silversmith W, Baker C, Deutsch D, Encarnacion-Rivera L, Kumar S, Burke A, Bland D, Gager J, Hebditch J, Koolman S, Moore M, Morejohn S, Silverman B, Willie K, Willie R, Yu SC, Murthy M, Seung HS. FlyWire: online community for whole-brain connectomics. Nat Methods. 2022 Jan;19(1):119-128. doi: 10.1038/s41592-021-01330-0. Epub 2021 Dec 23. PMID: 34949809; PMCID: PMC8903166.
- Schlegel P, Yin Y, Bates AS, Dorkenwald S, Eichler K, Brooks P, Han DS, Gkantia M, Dos Santos M, Munnelly EJ, Badalamente G. Whole-brain annotation and multi-connectome cell typing quantifies circuit stereotypy in Drosophila. BioRxiv. 2023 Jul 15.

- Matsliah A, et al. FlyWire Codex. Princeton Neuroscience Institute. https://codex.flywire.ai
- Cordella LP, et al. (2004). A (sub)graph isomorphism algorithm for matching large graphs. *IEEE TPAMI*, 26(10), 1367–1372.
