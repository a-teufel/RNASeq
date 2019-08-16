#! /bin/bash
#! This script runs a python script that adjusts the ratio of slow to fast codons from very small (approx. .05) to 1 and executes a new simulation for each ratio.

ratios=(1 10 25 50 75 90 100)
for i in "${ratios[@]}"
   do
     inner_dir="iterations/experiments/pp_rate/50_expression_error/weight_1/1st_iteration/"$i"ratio"
		 mkdir -p $inner_dir
		 cp -r 2_gene_extended_model/* $inner_dir
		 python 2_gene_extended_model/bin/edit_csv.py
		 pwd
		 echo $i | bc -l > $inner_dir/ratio.txt
   done

#Cleanup
cp 2_gene_extended_model/data/model_codons_reference.csv 2_gene_extended_model/data/model_codons.csv
