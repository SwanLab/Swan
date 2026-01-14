% LPV model
tora_lpv

% sampling intervals for each parameter
domain = [-45/180*pi 45/180*pi; -0.5 0.5];
% grid size: number of grid points for each parameter
gridsize = [23 23];

%% TP transformation, same as:
%   [S U] = tptrans(lpv, dep, domain, gridsize, 'close');

% sampling
lpvdata = sampling_lpv(lpv, dep, domain, gridsize);

% hosvd
[S U sv tol] = hosvd_lpv(lpvdata, dep, gridsize, 0.001);

% generating tight polytopic representation
hull = 'cno';
U = genhull(U, hull);
S = coretensor(U, lpvdata, dep);

% plot the results
plothull(U, domain);

% check model approximation error
[maxerr meanerr] = tperror(lpv, S, U, domain, 100);
disp('max and mean error:'); disp(maxerr); disp(meanerr);

save('tora_data', 'S', 'U', 'n', 'domain', 'gridsize');
