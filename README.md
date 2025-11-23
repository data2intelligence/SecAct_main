
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproducing the SecAct paper

These are scripts for creating the secreted protein signaling
signatures, validating their prediction accuracy, and demonstrating the
downstream applications of
<a href="https://github.com/data2intelligence/SecAct"
target="_blank">SecAct</a>.

## 1 Signature creation

Creation of secreted protein signatures.

## 1.1 Sample

    Rscript Secretome_s1_sample_summary.R

## 1.2 VST normalization

    swarm --verbose 1 Secretome_s2_coexpr_1_vst.swarm

## 1.3 Separate the signature

    swarm --verbose 1 Secretome_s2_coexpr_1_separate.swarm

## 1.4 Gating

    swarm --verbose 1 Secretome_s3_QC_1_gating.swarm

## 1.5 Summary

    Rscript Secretome_s3_QC_2_validate_summary.R

## 1.6 Combine the signature

    swarm --verbose 1 Secretome_s4_comb_0_all_in_one.swarm

## 1.7 Create composite signature

    swarm --verbose 1 Secretome_s4_comb_2_composite_signature.swarm
    Rscript Secretome_s4_comb_2_composite_signature_summary.R

## 1.8 Find the optial lambda

    swarm --verbose 1 Secretome_s4_model_1_explore_lambda.swarm
    Rscript Secretome_s4_model_2_explore_lambda_summary.R

## 2 Model validation

Validation using clinical data.

## 2.1 Anti-Netrin-1 therapy

    swarm --verbose 1 Secretome_s5_validation_1_NTN1_1_activity.swarm
    Rscript Secretome_s5_validation_1_NTN1_2_summary.R

## 2.2 Cytokine-blocking treatment

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

## 3 Downstream application

SecAct could be to three categories of input data: spatial
transcriptomics with multicellular or single-cell resolutions,
single-cell RNA-seq, and bulk transcriptomics.

