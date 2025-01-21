function y = multcoef(n,x)
% MULTCOEF Multinomial coefficient.
%  Y = MULTCOEF(N,X) returns the multinomial coefficient with parameter N at the values in X.
%  Let {X1, X2,  . . .  , Xk}, k > 1, be a set of random variables, each of which can take 
%  the values 0, 1,  . . .  , n; such that for every set of k nonnegative integers {n1,  
%  . . .  , nk} whose sum is n, the multinomial coefficient is, 
%
%                                                        (n1 + n2 + ... nk)!
%        C(n; n1, n2, ..., nk) = (n1, n2, ... nk)! =  -------------------------  .  
%                                                     n1! × n2! ×  . . .  × nk!    
%
%  It is possible to work with large factorials.
%
%  Syntax: y = multcoef(n,x) 
%      
%  Inputs:
%       n - number of trials.
%       x - vector of the interested values. 
%  Outputs:
%       y - multinomial coefficient.
%
%  Example. Assume that a die is thrown 60 times (n=60) and a record is kept of the number of 
%  times a 1, 2, 3, 4, 5, or 6 is observed. Each trial (e.g., throw of a die) has a (1 or 2 or
%  3 or . . . or 6) mutually exclusive outcome. In the die tossing data, k = 6, we are 
%  interested to get the multinomial coefficient.
%
%  Calling on Matlab the function: 
%             multcoef(n,x)
%
%  where n=60 and x=[13,10,8,10,12,7];
%
%  Answer is:
%
%  ans
%      = 1.0425e+042
%
%  Created by A. Trujillo-Ortiz, R. Hernandez-Walls and A. Castro-Perez
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%  Copyright (C) January 23, 2005.
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A., R. Hernandez-Walls and A. Castro-Perez. (2005). multcoef:
%    Multinomial coefficient. A MATLAB file. [WWW document]. URL http://
%    www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=6786
%
%  References:
% 
%  Abramowitz, M. and Stegun, I. A. (1964), Handbook of Mathematical
%           Functions, Government Printing Office, 823.24.1.2. Available on
%           Internet at the URL address http://hcohl.shell42.com/as/frameindex.htm
%

if nargin < 2, 
   error('You need to input two arguments.');
   return,
end;

if (length(n)~=1) | (fix(n) ~= n) | (n < 0),
   error('n must be a positive integer.');
   return,
end;

if sum(x) ~= n
    error('Inputs must satisfy n = x1 + x2 ... + xi.');
    return,
end;

factor1 = sum(log(1:n));

c = length(x);

f = [];
for i = 1:c,
    f = [f log(1:x(i))];
end;

factor2 = sum(f);

y = exp(factor1-factor2);
y = round(y);

return,
