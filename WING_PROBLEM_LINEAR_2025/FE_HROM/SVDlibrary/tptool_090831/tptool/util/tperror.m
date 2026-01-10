function [max_err, mean_err] = tperror(lpv, S, U, domain, n)
%TPERROR Calculate the error of the discretized TP model at random points
%	[max_err, mean_err] = TPERROR(lpv, S, U, domain, n)
%	
%	lpv      - LPV model
%	S        - core tensor of the TP model
%	U        - weight function data of the TP model
%	domain   - parameter domain of the LPV model
%	n        - number of random test points
%
%	max_err  - maximum L2 error over the n test points
%	mean_err - mean L2 error over the n test points

% TODO: use queryw instead of queryw1

% n random points in the given parameter space
P = size(domain,1);
x = zeros(n, P);
for i = 1:P
	x(:,i) = rand(n,1) * (domain(i,2) - domain(i,1)) + domain(i,1);
end

err = zeros(n, 1);
for i = 1:n
	% S from original model
	S_lpv = querylpv(lpv, x(i,:));
	
	% S from tp model
	W = queryw1(U, domain, x(i,:));
	S_tp = squeeze(tprod(S, W));
	
	% L2 error
	err(i) = norm(S_lpv - S_tp, 2);
	if err(i) > 100
		disp('Huge error')
		disp(x(i,:))
		S_lpv
		S_tp
	end
end

max_err = max(err);
mean_err = mean(err);
