function [min_kappa] = cal_min_kappa(theta,phifunct,g_nodal)

beta1 = sin((1-kappa)*theta)/sin(theta);
beta2 = sin(kappa*theta)/sin(theta);
phifunct = beta1*phifunct + beta2*g_nodal;

end

