function [N dN_dXi] = QuadElem_N_derN 

 
e = [-1 -1 1 1]; 
s = [-1 1 1 -1]; 
  % s = [-1 -1 1 1]; 
 % e = [-1 1 1 -1]; 


ngaus = 4; nstrain = 4; nnode = 4;  ndime = 2 ; 
N = zeros(ngaus,nnode);
dN_dXi = zeros(ndime*ngaus,nnode);
for g = 1:ngaus
    for a=1:nnode
        N(g,a) = 0.25*(1+s(a)*s(g)/sqrt(3))*(1+e(a)*e(g)/sqrt(3)); 
        dN_dXi(2*(g-1)+1:2*g,a) = 0.25*[(1+e(a)*e(g)/sqrt(3))*s(a)
                                        (1+s(a)*s(g)/sqrt(3))*e(a) ] ;
    end
end

 

 