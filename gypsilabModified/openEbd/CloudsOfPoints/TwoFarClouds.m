function[X,Y,V] = TwoFarClouds(M,N)
% [X,Y,V] = TwoFarClouds(M,N)
% inputs : M, N two integers
% outputs : 
% - X of size Mx2 points near the origin
% - Y of size Nx2 points near [10,10]
% - V real random vector of unit l^1 norm. 
X = randn(M,2);
Y = randn(N,2) + 10;
V = randn(N,1);
V = V/norm(V,1);
end