function [gfunc] =  cal_nodgfunc_objective_macro(igaus,work1,work2,phigp,element)


epsi  = element.material.opt_epsi;
%kappa = element.material.opt_kappa;
L     = element.material.opt_L;
mu = element.material.mu;
coeff = (lambda+3*mu)/(lambda+mu);

c1 = (epsi-1)/(coeff*epsi+1)*(coeff+1)/2;
c2 = (epsi-1)*(coeff-2)/(coeff + 2*epsi -1);
alpha1 = 2*c1;
alpha2 = c1*c2;

d1 = (1-epsi)/(coeff+epsi)*(coeff+1)/2;
d2 = (1-epsi)*(coeff-2)/(coeff*epsi + 2 -epsi);
beta1 = 2*d1;
beta2 = d1*d2;

g_menos = alpha1*work1 + alpha2*work2  ;
g_mas = -beta1*work1 - beta2*work2 ;


gfunc(igaus,:)= -g_menos;
outside= (phigp > 0);
if any(outside)
    gfunc(igaus,outside)= g_mas(outside);
end

end 