function [genefile, initfile, rnafile] = build_transcripts(filter_genes)
% [genefile, initfile, rnafile] = build_transcripts(filter_ids)
%   prepare fasta files for build_model(), including:
%   transcript sequences, initiation rates, mRNA levels.
%   currently using (Ingolia, 2009) ribo-seq and RNA-seq data.
%
% Input:
%   [filter_ids] (optional) gene indices / quantiles / names to include in
%   model. when quantiles (values in [0, 1]) are given, will select genes
%   based on RP data.
%
% Changelog:
%   08/03/17: the total mRNA copies of selected genes is normalized to
%             TOTAL_MRNA.
%
% February 2016
% Alon Diament / Tuller Lab

MEDIAN_INIT_RATE = 0.09;  % mRNA/s (Ciandrini, 2013)
TOTAL_MRNA = 30e3; % literature reports anything between 15,000 (von der Haar, 2008) and 60,000 (Zenklusen, 2008)
% MIN_MRNA_COUNT = 1;  % copies

load('data/sacCer', 'gene_id', 'gene_seq');
N = size(gene_id, 1);

if nargin < 1
    filter_genes = 1:size(gene_id, 1);
elseif iscell(filter_genes)
    % according to name, otherwise according to index / quantile
    filter_genes = any(ismember(gene_id, filter_genes), 2);
end

load('data/Ingolia_2009', 'RP', 'mRNA');

if any(filter_genes < 1 & filter_genes > 0)
    [~, ranks] = sort(RP);
    filter_genes = ranks(round(filter_genes * N));
end

init_rates = RP ./ mRNA;
init_rates(isnan(init_rates)) = 0;  % missing
init_rates = init_rates(filter_genes) / nanmedian(init_rates) * MEDIAN_INIT_RATE;

% mRNA_scale = min(mRNA(mRNA >0 )) / MIN_MRNA_COUNT;
mRNA_scale = sum(mRNA(filter_genes)) / TOTAL_MRNA;
mRNA = mRNA(filter_genes) / mRNA_scale;

gene_seq = gene_seq(filter_genes)';
gene_id = gene_id(filter_genes, 1);

genefile = 'data/model_genes.fasta';
initfile = 'data/model_init.csv';
rnafile = 'data/model_rna.csv';
%if exist(genefile, 'file')
    % avoid appending
%    delete(genefile);
%end

%fastawrite(genefile, struct('Header', gene_id, 'Sequence', gene_seq));
%writetable(table(gene_id, init_rates), initfile, 'Delimiter', '\t', ...
%           'WriteVariableNames', false);
%writetable(table(gene_id, mRNA), rnafile, 'Delimiter', '\t', ...
%           'WriteVariableNames', false);
