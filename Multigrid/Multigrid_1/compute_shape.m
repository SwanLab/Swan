% Store shape functions for rectangular elements
function [N] = compute_shape(zi,eta, nnode)
N = zeros(nnode,1);
N(1) = 0.25*(1-zi)*(1-eta);
N(2) = 0.25*(1+zi)*(1-eta);
N(3) = 0.25*(1+zi)*(1+eta);
N(4) = 0.25*(1-zi)*(1+eta);
end