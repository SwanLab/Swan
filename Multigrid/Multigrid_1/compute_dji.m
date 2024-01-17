function [dji] = compute_dji(xt,zi,eta, nnode, ndim)
temp_mat = compute_dj(xt,zi,eta, nnode, ndim);
num = 1/det(temp_mat);
dji = zeros(2,2);
dji(1,1) = temp_mat(2,2);
dji(2,2) = temp_mat(1,1);
dji(1,2) = -temp_mat(1,2);
dji(2,1) = -temp_mat(2,1);
dji = num*dji;
end