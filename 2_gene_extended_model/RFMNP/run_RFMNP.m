% toy example for the RFMNP, based on yeast genes
%
% NOTE: using fasta file for gene sequences, but it's possible to give
%       build_model() just a variable (table of [Genes]) and get [Elong].
%
% February 2016
% Alon Diament / Tuller Lab

chunksize = 10;  % codons
poolsize = 1200;  % ribosome pool

codonfile = build_codons();
[seqfile, initfile, rnafile] = build_transcripts([0.5, 0.75, 0.90, 0.95]);  % select transcripts according to quantiles
[Elong, Genes, Init, mRNA, errori] = build_model(codonfile, seqfile, ...
    chunksize, initfile, rnafile);
if errori
    error('build_model() failed (error %d)', errori);
end

ind = cellfun(@length, Elong) > 3 & mRNA > 0 & Init > 0;
Genes = Genes(ind);
Elong = Elong(ind);
Init = Init(ind);
mRNA = mRNA(ind);

[T, Rate, PA, Density, Pool] = solve_RFMNP(Elong, Init, mRNA, poolsize);
% figure;
% plot_RFMNP(T, Rate, Density, Pool, Init, mRNA)
