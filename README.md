
# FFPE-MicroDBiT

FFPE-MicroDBiT is a spatial technology for formalin-fixed paraffin-embedded (FFPE) tissues that enables the recovery of spatially resolved bacterial transcriptomes alongside host gene expression, revealing species-specific bacterial distributions and transcriptional heterogeneity.

## Overview

<img width="971" height="370" alt="FFPE-MicroDBiT" src="https://github.com/user-attachments/assets/95549361-c646-43eb-9faa-71344d7bf924" />


FFPE-MicroDBiT is a spatial transcriptomic platform that enables simultaneous profiling of microbial and host RNAs within the same FFPE tissue section. By integrating in situ RNA polyadenylation and random-primed transcript capture with deterministic microfluidic barcoding, followed by rRNA depletion, FFPE-MicroDBiT enables spatially resolved characterization of microbial transcriptional profiles alongside host gene expression in preserved tissues.

### Key Capabilities

- Taxonomic classification of host-depleted reads
- Spatially resolved bacterial species quantification
- Spatial characterization of host–microbial associations
- Application to archival human colon specimens

## Repository Structure

- **Data_processing/** — Scripts for sequencing data processing, taxonomic classification, and spatial expression matrix generation.
- **Data_analysis/** — Scripts for downstream computational and statistical analyses of host and microbial spatial transcriptomics data.
- **Visualization/** — Scripts for generating figures and visualizing spatial expression patterns and analysis results.

## Usage

Scripts can be executed independently or as part of the relevant analysis workflows. Users should update the input/output paths and parameters according to their datasets.

Required input data, reference files, and software dependencies vary among scripts.
