% Compute derivative of shape functions with respect to x and y directions
function [dfx] = compute_dfx(xt,zi,eta, nnode, ndim)
dfx = compute_dfxi(zi,eta, nnode, ndim)*compute_dji(xt,zi,eta, ...
    nnode, ndim);
end