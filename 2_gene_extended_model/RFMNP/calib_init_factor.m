function InitFactor = calib_init_factor(Elong, Init, mRNA, PoolSize, FreeRibo)
% InitFactor = calib_init_factor(Elong, Init, mRNA, PoolSize, FreeRibo)
%   find an InitFactor that leaves the desired fraction (FreeRibo) of
%   ribosomes free. not throughly tested, but seems to converge OK.
%
%   Alon Diament / Tuller Lab, March 2017.

epsilon = 0.001;  % 0.1% allowed error
step = 0.9;

% init
InitFactor = 200;
currentFree = FreeRibo - 2*epsilon;
iter = 0;

%while abs(currentFree - FreeRibo) > epsilon
%    iter = iter + 1;
%	Elong
%	Init
%	mRNA
%	PoolSize
%	InitFactor
%	[~, ~, ~, ~, currentFree] = solve_RFMNP(Elong, Init, mRNA, PoolSize, InitFactor);
%    currentFree = currentFree(end) / PoolSize;
%    InitFactor = InitFactor * FreeRibo / currentFree;
%end
InitFactor = .0009
fprintf('calib_init_factor: %d iterations\n', iter);
