function [h1,h2] = coef_DtC(gamma,element,problembsc)

switch problembsc.phisical_type
    case {'ELASTIC'}
        poiss_ref = element.material.poiss;
        young = element.material.young;
        alpha = (1+poiss_ref)./(1-poiss_ref);
        beta = (3-poiss_ref)./(1+poiss_ref);
        coeff1 = (1-gamma)./(1+alpha.*gamma);
        coeff2 = (1-gamma.*(alpha-2*beta))./(1+beta.*gamma);
        h1 = - 1./young.*coeff1*4;
        h2 = 1./young.*coeff1.*coeff2;
    case {'THERMAL'}
        k = element.material.k11;
        h1 = -k*(1-gamma)/(1+gamma);
        h2 = 0;
end

end