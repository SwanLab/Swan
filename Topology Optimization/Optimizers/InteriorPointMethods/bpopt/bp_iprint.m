function [] = bp_iprint(bp,iter,x,lam,zL,zU,alpha_pr,alpha_du,s,xL,xU,bL,bU)
n = size(x,2);
screen = 1;
if (mod(iter,15)==0),
   fprintf(screen,'   Iter       Merit   Objective   log10(mu)        Pcvg        Dcvg    alpha_pr    alpha_du\n');
end
du = sum(abs(bp_objgrad(bp,x,s)' + bp_jac(bp,x,bL,bU)'*lam' - zL' + zU'));
me = bp_merit(bp,x,xL,xU,s,bL,bU);
ob = bp_obj(bp,x);
fprintf(screen,'  %5i %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n',iter,me,ob,log10(bp.mu),bp_theta(bp,x,s,bL,bU),du,alpha_pr,alpha_du);
