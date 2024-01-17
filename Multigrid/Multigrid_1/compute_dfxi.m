% Store the derivative of shape function with respect to zi and eta
% directions, as a 2D array
function [dfxi] = compute_dfxi(zi,eta, nnode, ndim)
dfxi = zeros(nnode,ndim);
dfxi(1,1) = -0.25*(1-eta);
dfxi(2,1) = 0.25*(1-eta);
dfxi(3,1) = 0.25*(1+eta);
dfxi(4,1) = -0.25*(1+eta);
dfxi(1,2) = -0.25*(1-zi);
dfxi(2,2) = -0.25*(1+zi);
dfxi(3,2) = 0.25*(1+zi);
dfxi(4,2) = 0.25*(1-zi);
end