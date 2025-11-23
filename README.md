
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the SecAct paper

These are scripts for generating the signaling signatures of secreted
proteins, evaluating their prediction accuracy, and demonstrating the
downstream applications of
<a href="https://github.com/data2intelligence/SecAct"
target="_blank">SecAct</a>.

## 1 Signature creation

## 2 Model validation

Validation using in vivo data from secreted protein blockades

## 2.1 Anti-Netrin-1 therapy

    swarm --verbose 1 Secretome_s5_validation_1_NTN1_1_activity.swarm
    Rscript Secretome_s5_validation_1_NTN1_2_summary.R

## 2.2 Cytokine-blocking treaement

    swarm --verbose 1 Secretome_s5_validation_2_blocking_1_activity.swarm
    Rscript Secretome_s5_validation_2_blocking_2_summary.R

## 2.3 Anti-VEGF therapy

    swarm --verbose 1 Secretome_s5_validation_3_VEGF_1_activity.swarm
    Rscript Secretome_s5_validation_3_VEGF_2_summary.R

## 2.4 Pathway and proteomics in CPTAC patient cohorts

    swarm --verbose 1 Secretome_s5_validation_4_CPTAC_1_activity.R
    Rscript Secretome_s5_validation_4_CPTAC_2_summary.R

## 2.5 Pathway and proteomics in DepMap cell-line data

    swarm --verbose 1 Secretome_s5_validation_5_ProCan-DepMapSanger_1_activity.swarm
    Rscript Secretome_s5_validation_5_ProCan-DepMapSanger_2_summary.R

## 2.6 Transcription factor in scRNA-seq data

    swarm --verbose 1 Secretome_s5_validation_6_scRNAseq_0_preprocess.swarm
    swarm --verbose 1 Secretome_s5_validation_6_scRNAseq_1_activity.swarm
    Rscript Secretome_s5_validation_6_scRNAseq_2_summary.R

## 3. Downstream application

SecAct could be to three categories of input data: spatial
transcriptomics with multicellular or single-cell resolutions,
single-cell RNA-seq, and bulk transcriptomics.

