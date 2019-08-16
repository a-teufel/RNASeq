echo "input last iteration num"
read last_num
echo "input new iteration num"
read num
##replace old exp num with new
find bin/make_experiments.sh -type f -exec sed -i "s/"$last_num"_iteration/"$num"_iteration/g" {} \;
find 2_gene_extended_model/ProjectScript.m -type f -exec sed -i "s/"$last_num"_iteration/"$num"_iteration/g" {} \;
sed -e "s/$last_num/$num/g" bin/execute_$last_num\.sh > bin/execute_$num\.sh
chmod +x bin/execute_$num\.sh
