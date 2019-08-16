This experiment will be ran with the following save function every 1k generations:

save(name,'Elong','PoolR','Elong0','scoresList','allGenes','Density0','Density','Genes','Rate','Rate0','currentScore','k','mutatedScore','count','-v7')
csvwrite('../iterations/1st_iteration/csv_output/scores.csv',scoresList)
csvwrite('../iterations/1st_iteration/csv_output/current_elongation.csv',Elong)
csvwrite('../iterations/1st_iteration/csv_output/init_elongation.csv',Elong0)
csvwrite('../iterations/1st_iteration/csv_output/1_gene_codons.csv',Genes(1).Codons)
csvwrite('../iterations/1st_iteration/csv_output/2_gene_codons.csv',Genes(2).Codons)
csvwrite('../iterations/1st_iteration/csv_output/elongation_rates.csv','Elong')

2 Genes (large, small)

N = 100

W = 0.0

250K Generations
