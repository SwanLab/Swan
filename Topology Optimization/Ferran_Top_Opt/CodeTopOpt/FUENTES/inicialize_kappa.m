function kappa_ini = inicialize_kappa(design_variable,gradient,kappa_old,scalar_product,algorithm)

switch algorithm
    case 'level_set'
        kappa_ini = 1;
        
    case 'Projected_gradient'
        if isempty(kappa_old)
            norm_gamma = sqrt(scalar_product(design_variable,design_variable));
            norm_g = sqrt(scalar_product(gradient,gradient));
            kappa_ini = norm_gamma/norm_g;
        else
            kappa_ini = 10*kappa_old;
        end
end


end