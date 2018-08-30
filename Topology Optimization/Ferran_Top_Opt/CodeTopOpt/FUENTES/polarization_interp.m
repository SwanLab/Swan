function [P_fun,P] = polarization_interp(mu,K)
dmu = diff(mu);
%lam = 2*mu*(3*K-2*mu)/(3*K+4*mu); %lam
d = 2; % the dimension, 2 in 2D
lam = K - 2/d*mu;
dlam = diff(lam);

C = [2*mu+lam lam 0;lam 2*mu+lam 0;0 0 mu];
dC = [2*dmu+dlam dlam 0;dlam 2*dmu+dlam 0;0 0 dmu];

P = (C\dC);
P_fun = matlabFunction(P);
end