% LPV model
wing_lpv

filename='wing_data';

% intervals
domain = [15 45; -0.1 0.1];
% grid size
gridsize = [31, 31];


% TP transformation:
[S U] = tptrans(LPV, dep, domain, gridsize, 'cno');

% plot the results
plothull(U, domain);

% check model approximation error
[maxerr meanerr] = tperror(LPV, S, U, domain, 100);
disp('max and mean error:'); disp(maxerr); disp(meanerr);

save(filename, 'S', 'U', 'n', 'domain', 'gridsize');
