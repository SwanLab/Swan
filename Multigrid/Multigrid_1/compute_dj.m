function [dj] = compute_dj(xt, zi,eta, nnode, ndim)
dj = xt*compute_dfxi(zi, eta, nnode, ndim);
end