function dp_dt = dDensity_RFMNP(Density, Elong, Init, mRNA, InitFactor)
% dp_dt = dDensity_RFMNP(Density, Elong, Init, mRNA, C)
%   differential equations of the RFMNP (Raveh et al., 2016),
%   i.e. Ribosome Flow Model Network with Pool.
%   [Elong] is a cell array containing the rate vectors of all sites per gene.
%   [Init] are initiation rates.
%   densities are in [0, 1] for each site, but transcripts are weighted
%   according to given [mRNA] levels.
%   [InitFactor] is an initiation function (G) parameter.
%
% June 2015
% Alon Diament / Tuller Lab.
% (based on dDensity.m)
%
% Chagelog:
%   21/02/16: [Init]/[mRNA] transcript-specific weights for initiation
%             rates / mRNA levels.
%   08/03/17: [C] an initiation constant

ORF_size = cellfun(@length, Elong);
nORFs = length(ORF_size);
dp_dt = cell(nORFs, 1);

% ribosome pool
z = Density(1);  % free ribos pool
G = Init .* tanh(z/InitFactor);  % input per transcript
Density = mat2cell(Density(2:end), ORF_size + 1);  % 1 ORF in each cell
y = cellfun(@(x, y) x(end-1)*y(end), Density, Elong);  % outflowing ribos (see dp_dt(N+1))
dz_dt = sum(mRNA .* y) ...  % all outflowing ribos from transcripts
    - sum(mRNA .* G .* (1 - cellfun(@(x) x(1), Density)));  % all inflowing ribos (1st term in dp_pt(1))

for o = 1:nORFs
    % a single RFM
    N = ORF_size(o);
    dp_dt{o} = zeros(N+1, 1);
    dp_dt{o}(1) = G(o)*(1-Density{o}(1)) - Elong{o}(1)*Density{o}(1)*(1-Density{o}(2));  % init - out

    dp_dt{o}(2:(N-1)) = Elong{o}(1:(N-2)).*Density{o}(1:(N-2)).*(1-Density{o}(2:(N-1))) ...
        - Elong{o}(2:(N-1)).*Density{o}(2:(N-1)).*(1-Density{o}(3:N));  % in - out

    dp_dt{o}(N) = Elong{o}(N-1)*Density{o}(N-1)*(1-Density{o}(N)) ...
        - Elong{o}(N)*Density{o}(N);  % in - terminate
    dp_dt{o}(N+1) = Elong{o}(N)*Density{o}(N);  % this is the outflow y_j
end

% return a vector of derivatives
dp_dt = [dz_dt; cell2mat(dp_dt)];

% % test that no ribosome got lost
% ribo_test = mRNA .* cellfun(@(x) sum(x(1:end-1)), Density);  % all translating
% fprintf('%.2f, ', sum(ribo_test) + z);  % plus pool = all ribos (const)
