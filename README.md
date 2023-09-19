[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8357271.svg)](https://doi.org/10.5281/zenodo.8357271)

# Coupling of response biomarkers between tumour and peripheral blood in mesothelioma

### Introduction
This repository contains the computational workflows for the above-mentioned manuscript. Each folder represents a [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow. Broadly, individual steps are grouped by source data (bulk RNAseq, single cell RNAseq, Nanostring) in each folder. 

### Running the workflows

These workflows are containerised. Therefore, invoke the --run-singularity flag when running snakemake. The two sofware enviroments used across these workflows are large (5GB - 10GB). 

### Runtime statistics

The workflows were executed on a 16-core, 64GB RAM linux machine. Runtime statistics are included for bulk RNAseq, single cell RNAseq and Nanostring data as html reports auto-generated from snakemake.
