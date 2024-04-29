function [pospg,wpg] = Quadrature(ngaus) 
% [pospg,wpg] = Quadrature(ngaus) 
% pospg, wpg:   Gauss points and weights on the reference element
% ngaus:        number of Gauss points
%

if ngaus == 4 
   pos1 = 1/sqrt(3); 
   pospg = [-pos1   -pos1 
             pos1   -pos1 
             pos1    pos1 
            -pos1    pos1]; 
   wpg = [ 1 1 1 1]; 
else 
   error('Unavailable quadrature') 
end 
 
           
