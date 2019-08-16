echo "Enter last experiment weight: \n"
read old_weight
echo "Enter new experiment weight: \n"
read new_weight
find bin/make_experiments.sh -type f -exec sed -i "s/weight_"$old_weight"/weight_"$new_weight"/g" {} \;
find 2_gene_extended_model/ProjectScript.m -type f -exec sed -i "s/weight_"$old_weight"/weight_"$new_weight"/g" {} \;

bash -c "bin/make_experiments.sh"
bash -c "bin/make_iterations_from_config.sh"
