function PoolSize = calib_pool(Elong, Init, mRNA, InitFactor, FreeRibo)
% PoolSize = calib_pool(Elong, Init, mRNA, InitFactor, FreeRibo)
%   find a PoolSize that leaves the desired fraction (FreeRibo) of
%   ribosomes free. not throughly tested, but seems to converge OK.
%
%   Alon Diament / Tuller Lab, March 2017.

epsilon = 0.001;  % 0.1% allowed error
step = 0.9;

% init
PoolSize = sum(mRNA) / FreeRibo;
currentFree = FreeRibo - 2*epsilon;
iter = 0;

while abs(currentFree - FreeRibo) > epsilon
    iter = iter + 1;
    [~, ~, ~, ~, currentFree] = solve_RFMNP(Elong, Init, mRNA, PoolSize, InitFactor);
    currentFree = currentFree(end) / PoolSize;
    PoolSize = PoolSize + max(-step, step^iter * log2(FreeRibo / currentFree)) * PoolSize;
end

fprintf('calib_pool: %d iterations\n', iter);
