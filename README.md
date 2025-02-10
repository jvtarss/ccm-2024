# Microbial communities diversity analysis in microbial fuel cells

## Overview
Microbial fuel cells (MFCs) provide a sustainable solution for wastewater treatment and renewable energy generation, with their efficiency being highly dependent on microbial community composition, which varies considerably across effluent types and operational conditions. In this study, sequencing libraries of 16S rRNA amplicons were obtained from 227 samples across 31 published MFC studies.

## Key findings
- Dominant phyla: _Pseudomonadota_ (13.45%), _Bacteroidota_ (9.88%);
- Key electroactive genera: _Geobacter_ (1.10%), _Proteiniphilum_ (0.84%);
- Exclusive specialized metabolic pathways in different effluent conditions;
- Complex microbial interactions in bioelectrochemical systems.

## Limitations
- Fluctuations in sequencing depths among samples, leading to significant differences in feature coverage;
- Rarefaction removed large amounts of data but retained the maximum possible features for as many samples as possible;
- The analysis lacks sufficient sensitivity to correlate specific microbial communities with many operational conditions.

## Data processing
- Tools: QIIME2 2024.5.1, PICRUSt2 v2.5.3, SCNIC v0.6.6, biom-format v2.1.15, matplotlib v3.8.4, seaborn v0.12.2, re (built-in), pandas v2.2.2, numpy v1.26.4, font_manager v0.1.0, sklearn v1.4.2 
- Taxonomic classification with Greengenes2 2024.09 classifier

## Repository access
https://github.com/jvtarss/ccm-2024

