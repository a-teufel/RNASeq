#! /bin/bash
# Script to execute Freq Analysis

echo "Gene number: >"
read number
/usr/local/MATLAB/R2016a/bin/glnxa64/MATLAB -r "addpath(genpath('~/extended_model/2_gene_extended_model'));load('Variables9000.mat');multigen_speed_analysis(allGenes($number,:),$number,9000);"
