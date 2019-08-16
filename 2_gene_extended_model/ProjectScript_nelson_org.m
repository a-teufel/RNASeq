clearvars;
close all;
if ismac
			system('caffeinate -dims &'); % prevent computer from going into sleep
end
addpath(genpath('RFMNP'));
addpath(genpath('data'));
addpath(genpath('analysis'));
addpath(genpath('Results'));

%set model parameters
codonfile     = 'data/model_codons.csv';
ChunkSize     = 10;  % in codons
PoolSize      = 200e3;  % 187e3 +/- 56e3 (von der Haar, 2008)
FreeRibo      = 0.30; % demand pool with 15% free ribosomes based on endogenous genes (see below)
mut_per_gen   = 1;  % random mutations to introduce per generation
max_exp_error = 0.05; % 5% allowed error in expression in mutations
K_gen         = 2e4;  % generations
Fit           = .01;
count_update  = 0;
count_fail    = 0;
datetimestr   = datestr(datetime('now'));
date_format = 'mm_dd_yy';
%       ---------------------------------------------------------------------------------
%       ---------------------------------------------------------------------------------

%build/calibrate RFMNP
[seqfile, initfile, rnafile] = build_transcripts(linspace(95, 98, 6)/100);  % select genes & init transcripts once
[ElongR, Genes, Init, mRNA, errori] = build_model(codonfile, seqfile, ...
ChunkSize, initfile, rnafile);  % [Genes] is an array of RNASeq objects
if errori
		error('build_model() failed (error %d)', errori);
end

%assert(abs(sum(mRNA) - 30e3) < 1e-6);

InitFactor = calib_init_factor(ElongR, Init, mRNA, PoolSize, FreeRibo);

%%Calculate theoretical maximum rate
MaxRate = CalcMaxRate(codonfile, Genes, ChunkSize, Init, mRNA, PoolSize, InitFactor)

[~, RateR, PAR, DensityR, PoolR] = solve_RFMNP(ElongR, Init, mRNA, PoolSize, InitFactor);
[Genes, realScore] = ScoreFn(Genes, RateR, MaxRate);
%       ---------------------------------------------------------------------------------
%       ---------------------------------------------------------------------------------

%create unique save location and copy README.txt, create subfolders in the experiment variables directory for our later analysis
%disp(pwd)
fileID = fopen('ratio.txt','r');
formatSpec = '%d %f';
sizeA = [1];
A = fscanf(fileID,formatSpec,sizeA)
fclose(fileID);
folder_name = strcat('/stor/work/Wilke/ateufel/RNASeq/RNASeq/2_gene_extended_model/1st_iteration/',num2str(A), '_ratio_',datestr(datetime('now'),date_format));
csv_folder = strcat('/stor/work/Wilke/ateufel/RNASeq/RNASeq/iterations/data/pp_rate/50_expression_error/weight_1/1st_iteration/',num2str(A), 'ratio/');
mkdir(folder_name);
mkdir(csv_folder);
copyfile('README.md',folder_name);
ratios_folder = strcat(folder_name,'/Plots/Ratios');
speed_folder  = strcat(folder_name,'/Plots/Speed');
mkdir(speed_folder);
mkdir(ratios_folder);
%       ---------------------------------------------------------------------------------
%       ---------------------------------------------------------------------------------

%start experiment
for i=1:1 %Total Experiments
		%initiate variables
		scoreDiffVector = [];
		updatedGenerations = [];
		allGenes = [];
		fitGenes = [];
		scoresList = [];
		fitGenerations = [];
		count = 0;
		% init sequence and scores
		nG = length(Genes);
		for g = 1:nG
				rng('shuffle')
				Genes(g) = RNASeq(aa2nt(Genes(g).AA), Genes(g).Name);  % Creates a random AA-preserving sequence with uniform codon bias
		end
		Elong0 = build_model(codonfile, Genes, ChunkSize);  % update model parameters according to sequence
		[~, Rate0, PA0, Density0, Pool0] = solve_RFMNP(Elong0, Init, mRNA, PoolSize,InitFactor);
		[Genes, currentScore] = ScoreFn(Genes, Rate0, MaxRate);  % per gene & total score
		initScore = currentScore;
    [ramp_1,bod_1] = first_codon_fraction(Genes(1),{'TAT','TAC'});
    [ramp_2,bod_2] = first_codon_fraction(Genes(2),{'TAT','TAC'});
		%       ---------------------------------------------------------------------------------
		%       ---------------------------------------------------------------------------------
		%save initial values
		name = strcat(folder_name,'/Variables0');
		save(name);
		fid = fopen(strcat(csv_folder,'1_gene_elong.dat'), 'w');
		fprintf(fid, '%i,', Elong0{1});
		fprintf(fid,'\n');
		fclose(fid);
		fid = fopen(strcat(csv_folder,'2_gene_elong.dat'), 'w');
		fprintf(fid, '%i,', Elong0{2});
		fprintf(fid,'\n');
		fclose(fid);
		fid = fopen(strcat(csv_folder,'1_gene_codons.dat'), 'w');
		fprintf(fid, '%s,', Genes(1).Codons{:});
		fprintf(fid,'\n');
		fclose(fid);
		fid = fopen(strcat(csv_folder,'2_gene_codons.dat'), 'w');
		fprintf(fid, '%s,', Genes(2).Codons{:});
		fprintf(fid,'\n');
		fclose(fid);
		fid = fopen(strcat(csv_folder,'1_gene_density.dat'), 'w');
		fprintf(fid, '%i,', Density0{1});
		fprintf(fid,'\n');
		fclose(fid);
		fid = fopen(strcat(csv_folder,'2_gene_density.dat'), 'w');
		fprintf(fid, '%i,', Density0{2});
		fprintf(fid,'\n');
		fclose(fid);
		fid = fopen(strcat(csv_folder,'pool.dat'), 'w');
		fprintf(fid, '%i,', Pool0);
		fprintf(fid,'\n');
		fclose(fid);
		fid = fopen(strcat(csv_folder,'1_gene_codon_ratio.dat'), 'w');
		fprintf(fid, 'Ramp: %f, Total: %f',ramp_1,bod_1);
		fprintf(fid,'\n');
		fclose(fid);
		fid = fopen(strcat(csv_folder,'2_gene_codon_ratio.dat'), 'w');
		fprintf(fid, 'Ramp: %f, Total: %f',ramp_2,bod_2);
		fprintf(fid,'\n');
		fclose(fid);

		for k = 1:K_gen  %Total Generations
				%Save variables and print progress
				if ~mod(k, 100)
					%Recalculate the ratio of fast to slow codons
						[ramp_1,bod_1] = first_codon_fraction(Genes(1),{'TAT','TAC'});
						[ramp_2,bod_2] = first_codon_fraction(Genes(2),{'TAT','TAC'});
						count_fail = 0;
						count_update = 0;
						name = strcat(folder_name,'/Variables',num2str(k));
						for i=1:length(allGenes)
								allGenes(i).Bases = '';
						end
						fid = fopen(strcat(csv_folder,'1_gene_elong.dat'), 'a');
						fprintf(fid, '%i,', Elong{1});
						fprintf(fid,'\n');
						fclose(fid);
						fid = fopen(strcat(csv_folder,'2_gene_elong.dat'), 'a');
						fprintf(fid, '%i,', Elong{2});
						fprintf(fid,'\n');
						fclose(fid);
						fid = fopen(strcat(csv_folder,'1_gene_codons.dat'), 'a');
						fprintf(fid, '%s,', Genes(1).Codons{:});
						fprintf(fid,'\n');
						fclose(fid);
						fid = fopen(strcat(csv_folder,'2_gene_codons.dat'), 'a');
						fprintf(fid, '%s,', Genes(2).Codons{:});
						fprintf(fid,'\n');
						fclose(fid);
						fid = fopen(strcat(csv_folder,'1_gene_density.dat'), 'a');
						fprintf(fid, '%i,', Density{1});
						fprintf(fid,'\n');
						fclose(fid);
						fid = fopen(strcat(csv_folder,'2_gene_density.dat'), 'a');
						fprintf(fid, '%i,', Density{2});
						fprintf(fid,'\n');
						fclose(fid);
						fid = fopen(strcat(csv_folder,'pool.dat'), 'a');
						fprintf(fid, '%i,', PoolR);
						fprintf(fid,'\n');
						fclose(fid);
						fid = fopen(strcat(csv_folder,'1_gene_codon_ratio.dat'), 'a');
						fprintf(fid, 'Ramp: %i, Total: %f',ramp_1,bod_1);
						fprintf(fid,'\n');
						fclose(fid);
						fid = fopen(strcat(csv_folder,'2_gene_codon_ratio.dat'), 'a');
						fprintf(fid, 'Ramp: %i, Total: %f',ramp_2,bod_2);
						fprintf(fid,'\n');
						fclose(fid);
						csvwrite(strcat(csv_folder,'scores_list.dat'),scoresList);
						csvwrite(strcat(csv_folder,'current_gen.dat'),k);
						csvwrite(strcat(csv_folder,'accepted_mutations.dat'),count);
				end
				%       ---------------------------------------------------------------------------------
				%       ---------------------------------------------------------------------------------
				%Mutates
				mutatedSeq = copy(Genes);  % init
				g = randi(nG, 1, 1)  % random mutations per generation
				mutatedSeq(g) = mutate(Genes(g));  % Creates a new, 'mutated' sequence by moving a random codon of the original MRNA sequence to a random index

				%Plugs into RFMNP Model
				Elong = build_model(codonfile, mutatedSeq, ChunkSize);
				[~, Rate, PA, Density, Pool] = solve_RFMNP(Elong, Init, mRNA, PoolSize,InitFactor);
				%Update the new RNASeq with it's new measurements
				for n = 1:length(mutatedSeq)
						mutatedSeq(n).Elong    = Elong{n};
						mutatedSeq(n).Pool     = Pool;
						mutatedSeq(n).Density  = Density{n}(:,end);
						mutatedSeq(n).Rate     = Rate{n};
				end

				[mutatedSeq, mutatedScore] = ScoreFn(mutatedSeq, Rate, MaxRate);

				% test mutations further
				%random = .15*normrnd(0,.33)+.7;
				random = rand(1);
				Fit = RNASeq.FitFn(currentScore, mutatedScore,mRNA);
				if(random < Fit)
						% accept all mutations
						count=count+1;
						count_update = count_update + 1;
						Genes = copy(mutatedSeq);
						currentScore = mutatedScore;
				end
				%save the currentScore and gene for every generation into one vector
				allGenes = [allGenes, Genes];
				scoresList = [scoresList currentScore];
				end

				%       ---------------------------------------------------------------------------------
				%       ---------------------------------------------------------------------------------
				%Prints final results
				fprintf('exp %d: \tinit_score=%.2f \tfinal_score=%.2f \tfinal_fit=%.2f\n', ...
				i, initScore, currentScore, Fit);
				fprintf('\tinit_rate=%.2f \t\tfinal_rate=%.2f\n', ...
				mean(cellfun(@(x, y) x(end)/y(end), Rate0, RateR)), ...
				mean(cellfun(@(x, y) x(end)/y(end), Rate, RateR)));  % relative to real genes
				fprintf('\tinit_pool=%.2f \t\tfinal_pool=%.2f\n', ...
				Pool0(end)/PoolR(end), Pool(end)/PoolR(end));  % relative to real genes
				FreqAnalysis(fitGenes,'fitGenes');
				plot_ramp(Elong, Elong0, 'elongation rates', 'init');  % vs. init
				plot_ramp(Density, Density0, 'densities', 'init');
				pause(0.5);
end
