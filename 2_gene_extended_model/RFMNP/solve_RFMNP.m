function [T, Rate, PA, Density, Pool] = solve_RFMNP(Elong, Init, mRNA, Pool, InitFactor)
% [T, Rate, PA, Density, Pool] = solve_RFMNP(Elong, Init, mRNA, Pool)
%   This function solves the differential equation for RFMNP (Raveh et al.,
%   2016)and returns the protein production [Rate], the total amount of
%   protein produced [PA], the ribosome [Density] vector for all genes and
%   the free ribos in the [Pool] at times [T].
%   [Rate{i}(end)], [PA{i}(end)] and [Density{i}(1:N,end)] are the steady
%   state values.
%
% Input:
%   [Elong] internal elongation rates vector.
%   [Init] initiation rates vector.
%   [mRNA] mRNA abundance vector.
%   [Pool] pool of free ribosomes.
%   [InitFactor] an initiation function (G) factor.
%   (transcripts are empty from ribos at T_0)
%
% Hadas Zur, Alon Diament / Tuller Lab
% (adapted from dLam.m)
%
% Chagelog:
%    04/06/15: multiple RNAs with a finite ribosome pool (RFMNP). (Alon)
%    21/02/16: [Init]/[mRNA] transcript-specific weights for initiation /
%              mRNA level. (Alon)

ORF_size = cellfun(@length, Elong);
N = sum(ORF_size + 1);

% we need to start from a time interval T_int in which we
% solve our deferential equation. We make a guess regarding the interval
% (settling into equilibrium)

T_int = [0, max(ORF_size)*max(1./cellfun(@min, Elong))];

% We start from an empty ststem.

initial = zeros(N+1,1);
initial(1) = Pool; % ribosome pool

% We let Matlab solve our differential equation

% fprintf('pool test: ');
sol = ode15s(@(t, Density) dDensity_RFMNP(Density, Elong, Init, mRNA, InitFactor), ...
        T_int, initial);
% fprintf('\n');

% We return values along time

T=sol.x;
Density=sol.y;
Pool = Density(1, :);
Density = mat2cell(Density(2:end, :), ORF_size + 1); % 1 ORF in each cell
L = length(ORF_size);
Rate = cell(L, 1);
PA = cell(L, 1);

%L= Number of genes = length(cellfun(@length, Elong))
for o = 1:L
    Rate{o} = mRNA(o) * Density{o}(end-1, :) * Elong{o}(end);
    PA{o} = mRNA(o) * Density{o}(end, :);
end

end
