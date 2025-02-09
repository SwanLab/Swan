function [Nmatrix] = multinomial_powers_recursive(pow,ndim)

% computes the multinomial expansion of 
% (x_0 + x_1 + x_2 + ... + x_ndim)^pow 
% Nmatrix is a matrix of powers
% This is equivalent to finding all multi-indices with norm=1

% Need another thing to calculate the coefficients, but that is easy

% recursive on dimension!

if ndim==1,
  Nmatrix = pow;
else
  % recurse
  Nmatrix = [];
  for pow_on_x1 = 0:pow,
    % say we fix the power in the first dimension to be "pow_on_x1" (0,1,2,...)
    % then the possible terms are all terms for [(pow-pow_on_x1),ndim-1]
    [newsubterms] = multinomial_powers_recursive(pow-pow_on_x1,ndim-1);
    % stick on the power for the x1 part and add to Nmatrix
    Nmatrix = [Nmatrix; [pow_on_x1*ones(size(newsubterms,1),1) , newsubterms] ]; 
  end
end

