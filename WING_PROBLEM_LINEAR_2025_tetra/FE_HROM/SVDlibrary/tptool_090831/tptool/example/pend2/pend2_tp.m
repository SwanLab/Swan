pend2_lpv

% parameters: a_1, a_2, da_1, da_2
domain = [-pi/12 pi/12; -pi/12 pi/12; -pi pi; -pi pi];
gridsize = [81 81 21 21];

%[S U] = tptrans(LPV, dep, domain, gridsize, 'close');

reply = input('Sampling and HOSVD? Y/N [Y]: ', 's');
if isempty(reply) || lower(reply)=='y'
	% sampling
	lpvdata = sampling_lpv(LPV, dep, domain, gridsize);

	% hosvd
	[Scan Ucan sv tol] = hosvd_lpv(lpvdata, dep, gridsize, 0, [5 5 2 2]);

	save('pend2_data', 'lpvdata', 'dep', 'Scan', 'Ucan', 'n', 'domain', 'gridsize');
else
	load('pend2_data', 'lpvdata', 'dep', 'Scan', 'Ucan', 'n', 'domain', 'gridsize');
end

reply = input('Hull generation? Y/N [Y]: ', 's');
if isempty(reply) || lower(reply)=='y'
	% generating tight polytopic representation
	U = genhull(Ucan, 'cno');
	S = coretensor(U, lpvdata, dep);

	save('pend2_data', '-append', 'S', 'U');
else
	load('pend2_data', 'S', 'U');
end

% plot the results
plothull(U, domain);

% check model approximation error
[maxerr meanerr] = tperror(LPV, S, U, domain, 100);
disp('max and mean error:'); disp(maxerr); disp(meanerr);
