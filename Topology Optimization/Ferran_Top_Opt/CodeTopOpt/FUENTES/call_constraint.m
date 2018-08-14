function [constr,dconstr_gp] = call_constraint(gamma,dim,element,problembsc,coordinates,Msmooth)

Vfrac = element.material.Vfrac;

[vol,dvol_n(:,1)] = cal_omega(gamma,dim,element,problembsc,coordinates);
        
constr = ((vol/Vfrac - 1));
dconstr_gp = dvol_n/Vfrac;


end