#!/bin/bash

module load R/4.3

Rscript Secretome_s3_QC_1_gating.R $1 $2
