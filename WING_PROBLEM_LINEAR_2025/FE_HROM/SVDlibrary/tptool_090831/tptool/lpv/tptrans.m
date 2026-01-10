function [S U] = tptrans(lpv, dep, domain, gridsize, hull)
%TPTRANS TP model transformation
%	[S U] = TPTRANS(lpv, dep, domain, gridsize, hull)
%
%	lpv      - LPV model (cell array of functions)
%	dep      - parameter dependencies for each element of the LPV model
%	domain   - parameter domains
%	gridsize - sampling gridsize for each parameter
%	hull     - type of polytopic representation
%
%	S        - core tensor
%	U        - discretized weighting functions
%
%	See also SAMPLING_LPV, HOSVD_LPV, GENHULL, CORETENSOR


% TODO: matlab docs

% sampling
data = sampling_lpv(lpv, dep, domain, gridsize);
% hosvd
[Scan Ucan] = hosvd_lpv(data, dep, gridsize, 1e-6);
% generating polytopic representation
U = genhull(Ucan, hull);
S = coretensor(U, data, dep);
%for i = 1:length(U)
%	Ucan{i} = pinv(U{i})*Ucan{i};
%end
%S = tprod(Scan, Ucan);

%plothull(U, domain);
%tperror(lpv, S, U, domain, 1000)
