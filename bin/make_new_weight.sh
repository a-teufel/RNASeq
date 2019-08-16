echo "Enter last experiment weight: \n"
read old_weight
echo "Enter new experiment weight: \n"
read new_weight
mkdir -p iterations/experiments/pp_rate/50_expression_error/weight_$new_weight

#Edit the control directory
find 2_gene_extended_model/ProjectScript.m -type f -exec sed -i "s/weight_"$old_weight"/weight_"$new_weight"/g" {} \;
find 2_gene_extended_model/RNASeq.m -type f -exec sed -i "s/W = "$old_weight"/W = "$new_weight"/g" {} \;

#Edit the build scripts
find bin/make_experiments.sh -type f -exec sed -i "s/weight_"$old_weight"/weight_"$new_weight"/g" {} \;
find bin/execute_1st.sh -type f -exec sed -i "s/weight_"$old_weight"/weight_"$new_weight"/g" {} \;

bash -c "bin/make_experiments.sh"
echo "Made first iteration"
pwd
bash -c "bin/make_iterations_from_config.sh experiment_iterations.cfg"
