function [Nmatrix,Ncoef] = multinomial_expand(pow,ndim)
% MULTINOMIAL_EXPAND determines the matrix of powers for a multinomial expansion
% of the form (x_1 + x_2 + x_3 + ... + x_ndim)^pow
% 
% Nmatrix - matrix of powers, each row representing a single term in the expansion
%  for example, the row [0,1,0,2] would represent (x_2)*(x_4)^2
%  Note, this is equivalent to finding all multiindices 
%  k = [k_1,k_2,...,k_ndim] with |k|=sum(k)=pow
% 
% Ncoef - vector of coefficients (multinomial coefficient)
%
% Ex. Evaluating the multinomial at a point x = [x_1, x_2, x_3, ... ]
% -> sum(Ncoef .* repmat(x,size(Nmatrix,1),1).^Nmatrix)
% 
% Ex. Compute all multiindices of length 6 and order 4:
% -> Nmatrix = multinomial_expand(4,6); 
%
% 
% The method is recursive, so it can be a bit slow, but
% it can work with relatively large size inputs (e.g. ndim = 50,100,...)
% compared to previous versions.  This is because we only compute the 
% needed number of terms and don't use kron or factorial
%
% By Isaac Asher
% Inspired by multiNomial.m by Alabi Bojesomo
% and multinomial.m by Mukhtar Ullah

% recursive function for doing the matrix of powers
% equivalently, this computes all multi-indices k with |k|=pow
% (here k is the multiindex (k_1,k_2,...) and |k| = sum(k_i))
[Nmatrix] = multinomial_powers_recursive(pow,ndim);

% do coefs separately
if nargout>1,
  % inspired by multinomial.m by Mukhtar Ullah on the Matlab FileExchange
  powvec = repmat(pow,size(Nmatrix,1),1); % vector of powers, "n" in multinomial.m
  % Nmatrix is "k" in multinomial.m
  Ncoef = floor(exp(gammaln(powvec+1) - sum(gammaln(Nmatrix+1),2))+0.5);
end

