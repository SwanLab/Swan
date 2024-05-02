function [pospg, wpg] = Quadrature_cont(ngaus1d, side);
% [pospg, wpg] = Quadrature_cont(ngaus1d, side);
% Gauss points and weigths for the numerical quadrature on a side of an element
%

if ngaus1d == 2
    pos1 = 1/sqrt(3); 
    if side == [1,2]
        pospg = [-pos1   -1
                  pos1   -1];
    elseif side == [2,3]
        pospg = [1   -pos1
                 1    pos1];             
    elseif side == [3,4]
        pospg = [ pos1   1
                 -pos1   1];
    elseif side == [4,1]
        pospg = [-1    pos1
                 -1   -pos1];
    end
    wpg=[ 1 1 ]; 
else
    error('Unavailable quadrature on the boundary');   
end
