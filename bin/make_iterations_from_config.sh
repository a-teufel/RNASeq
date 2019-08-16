experiments=()
while IFS='' read -r line || [[ -n "$line" ]]; do
	experiments+=("$line")
done < "$1"
echo ${#experiments[@]}-1 >&2
for ((experiment=0; experiment< ${#experiments[@]}-1; experiment++ )); do
	if [ $experiment == ${#experiments[@]} ]; then
		echo $experiment'last one' >&2
		break
	fi
	echo $experiment >&2
	last_num=${experiments[experiment]}
	first_num=${experiments[0]}
	num=${experiments[experiment+1]}
	echo $num >&2
	##replace old exp num with new
	find bin/make_experiments.sh -type f -exec sed -i "s/"$last_num"_iteration/"$num"_iteration/g" {} \;
	find 2_gene_extended_model/ProjectScript.m -type f -exec sed -i "s/"$last_num"_iteration/"$num"_iteration/g" {} \;
	bash -c "bin/make_experiments.sh"
	sed -e "s/$last_num/$num/g" bin/execute_$last_num\.sh > bin/execute_$num\.sh
	chmod +x bin/execute_$num\.sh
done
#Cleanup
find bin/make_experiments.sh -type f -exec sed -i "s/"$num"_iteration/"$first_num"_iteration/g" {} \;
find 2_gene_extended_model/ProjectScript.m -type f -exec sed -i "s/"$num"_iteration/"$first_num"_iteration/g" {} \;
