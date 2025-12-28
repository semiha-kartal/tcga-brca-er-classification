# tcga-brca-er-classification

## Overview
This project implements a supervised machine learning pipeline to classify estrogen receptor (ER) status in breast cancer using RNA-seq gene expression data from the TCGA-BRCA cohort.

The aim of the project was to gain hands-on experience with end-to-end machine learning workflows for high-dimensional transcriptomic data, including data preprocessing, feature selection, model training, and performance evaluation.

## Data
- Gene expression: TCGA-BRCA RNA-seq (FPKM-UQ)
- Clinical labels: ER status by immunohistochemistry
- Source: TCGA via UCSC Xena / GDC

Only primary tumor samples were included.

## Methods
- Sample harmonization between expression and clinical data
- Removal of duplicated genes (keeping highest expressed transcript)
- Variance-based feature selection (top 3000 genes)
- Explicit prevention of label leakage (ESR1 excluded)
- Elastic net logistic regression (glmnet)
- Cross-validation and ROC/AUC-based evaluation

## Results
The model achieved robust classification performance on held-out test data.
The elastic net approach selected a sparse set of predictive genes, enabling biological interpretability.

## Motivation
This project was carried out independently to strengthen practical experience with supervised machine learning on real cancer transcriptomic data, bridging computational methods and cancer biology.

## Requirements
- R (>= 4.2)
- glmnet
- pROC
