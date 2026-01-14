function [z,w,Lambda,DATAOUT] = GET_ecmPOINTS(Xf,W,DATA)
%
% Empirical cubature method
% -----------------------------
% Given Xf and W, GET_ecmPOINTS returns a set of indices z
% and associated positive weights w so that the integration error of each
% column of Xf (as well as the integral of the domain)  is minimized: 
%     min(Xf'*W - Xf(z,:)'w  + (sum(W) -sum(w)))
%
% INPUTS
% ------
% Xf: M x P n matrix containing, columnwise, the integrand defined at the M
% integration points of the underlying grid (n is the number of entries of the  integrand, and P the number of
% samples)
% W: M x 1 vector containing the integration weights associated with the
% underlying grid 
%  Optional inputs 
% -------------------
%DATA.TOLsvdXf =0; % Truncation tolerance determining dominant modes of Xf (in %)
% of Xf . Default value = 1e-12
% 
% Joaquín A. Hernández, January 8-th - 2016/Nov-2021
% See paper: Dimensional hyperreduction of nonlinear parameterized finite element models using optimized cubature
%  hernandez2016dimensional.pdf
%---------------------------------------
%
if nargin == 0
 
    load('tmp2.mat')
        
end
