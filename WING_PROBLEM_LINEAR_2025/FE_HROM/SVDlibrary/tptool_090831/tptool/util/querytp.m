function Sp = querytp(S, U, domain, p)
%QUERYTP Query TP model at a given parameter point
%	Sp = QUERYTP(S, U, domain, p)
%	
%	S        - core tensor of the TP model
%	U        - (discrete) weight function data of the TP model
%	domain   - intervals for each dimension
%	p        - parameter point where we query the model
%
%	Sp       - system at p (linear interpollated weights)

% TODO: method param to interp1 (spline,..)

W = queryw1(U, domain, p);
Sp = squeeze(tprod(S, W));
