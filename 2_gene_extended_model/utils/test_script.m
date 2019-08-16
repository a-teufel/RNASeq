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
K_gen         = 2.5e4;  % generations
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

[~, RateR, ~, DensityR, PoolR] = solve_RFMNP(ElongR, Init, mRNA, PoolSize, InitFactor);
[Genes, realScore] = ScoreFn(Genes, RateR, PoolR, RateR);
array_rand = [0,1,1,1,1,1,2,3,4,5,6,1]
RNASeq.write_csv(Genes(1),'RNASeq1.csv',array_rand)
