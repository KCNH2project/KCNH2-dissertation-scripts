# KCNH2 Variant Interpretation Pipeline

This repository contains the custom R scripts used for data processing, annotation, evidence integration, and visualisation during an MSc Genomic Medicine research project investigating automated ACMG/ACGS re-classification of KCNH2 missense variants associated with Long QT Syndrome type 2 (LQT2).

The pipeline integrates multiple evidence streams relevant to clinical variant interpretation, including:

population frequency assessment using gnomAD filtering allele frequencies (FAF)
computational prediction evidence (REVEL and SpliceAI)
mutational hotspot mapping (PM1)
automated patch-clamp (APC) functional evidence integration
Oxford laboratory classification extraction and review
ACMG/ACGS evidence collation and Bayesian point-based re-classification

# Repository Structure

The scripts are designed to run sequentially as a modular analysis pipeline:

Script	Purpose
SCRIPT-1	Multi-laboratory variant dataset collation and cleaning
SCRIPT-2	HGVS nomenclature harmonisation and Ensembl Variant Recoder querying
SCRIPT-3	APC dataset import and processing
SCRIPT-4	gnomAD population frequency and FAF evidence assessment
SCRIPT-5	Oxford Alamut classification extraction
SCRIPT-6	Computational prediction score retrieval (REVEL, SpliceAI, VEP)
SCRIPT-7	PM1 mutational hotspot assessment
SCRIPT-8	APC functional evidence integration
SCRIPT-9	ACMG/ACGS evidence collation and final re-classification
SCRIPT-10	Figure generation and visualisation outputs

# Reference Sequences

Variant nomenclature harmonisation and annotation were performed against the canonical KCNH2 RefSeq transcripts:

NM_000238.4
NP_000229.1

# Software Environment

Analyses were performed primarily in R using custom scripts and external annotation resources including:

Ensembl Variant Effect Predictor (VEP)
gnomAD
REVEL
SpliceAI
Alamut
CardioDB tools

Some annotation steps requiring VEP plugins and associated resources were executed on an in-house laboratory R server environment.

# Notes

This repository contains scripts only and does not include identifiable patient data or restricted clinical datasets.

Certain external datasets, laboratory resources, and proprietary software tools used during the project (e.g. Alamut database access) are not included within this repository.
