function codonfile = build_codons()
% build_codons()
%   computing codon decoding times based on tAI weights.
%
% February 2016
% Alon Diament / Tuller Lab

MEDIAN_CODON_RATE = 10;  % aa/s (Arava, 2003)

load('data/tGCN.mat');  % tRNA gene copy numbers
stop_list = {'TGA', 'TAA', 'TAG'};

[W, codon] = calc_tAI_weights(tGCN);
W(ismember(codon, stop_list)) = [];
codon(ismember(codon, stop_list)) = [];

W = 1./W;  % times
W(25) = max(W(W~=W(25))) * 10;  % CGA fix (very slow)
W = W / median(W) / MEDIAN_CODON_RATE;

codonfile = 'data/model_codons.csv';
%writetable(table(codon, W), codonfile, 'Delimiter', '\t', ...
%           'WriteVariableNames', false);


function [W, codon_list] = calc_tAI_weights(tGCN, S, S_rules)
% [tGCN]: tRNA Gene Copy Numbers given according to anti-codons.
% [S] s-values for all possible interactions between codons and tRNAs
% [S_rules] all interaction rules, 1st anti-codon position (1st col) and 
%   3rd codon position (2nd col).

if nargin < 3
    S_rules = {'A',    'T'; ... % watson-crick
               'G',    'C'; ...
               'T',    'A'; ...
               'C',    'G'; ...
               'G',    'T'; ... % wobble
               'A',    'C'; ...
               'A',    'A'; ...
               'T',    'G'};
end
if nargin < 2
    S = [0, 0, 0, 0, 0.561, 0.28, 0.9999, 0.68]; % 0.41 original / 0.561 optimization vs. PA
end

t_GCN = tGCN.GCN;
t_anti = tGCN.anti_codon;
t_codon = cellfun(@seqrcomplement, t_anti, ...
    'UniformOutput', false);

codon_list = fieldnames(codoncount(''));
stop_list = {'TGA', 'TAA', 'TAG'};
nC = length(codon_list);
W = zeros(nC, 1);

for c = 1:nC
    if any(strcmp(codon_list{c}, stop_list))
        W(c) = NaN;
        continue;
    end
    % candidate anti-codons (compare first 2 pos in codon)
    cand = find(strncmp(codon_list{c}, t_codon, 2))';
    if isempty(cand)
        continue;
    end
    % check S-rules for 3rd position
    for cand = cand
        is_tRNA_fits = @(x, y) (t_anti{cand}(1) == x) & ... % 1st anti pos
            (codon_list{c}(3) == y); % 3rd codon pos
        idx = cellfun(is_tRNA_fits, ... 
            S_rules(:, 1), S_rules(:, 2));
        if any(idx)
            W(c) = W(c) + (1 - S(idx))*t_GCN(cand);
        end
    end
end

W(strcmp('ATG', codon_list)) = t_GCN(strcmp('ATG', t_codon));  % known
W = W / max(W);
W(W==0) = geomean(W(W~=0 & ~isnan(W)));
