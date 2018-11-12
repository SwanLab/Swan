function design_variable_new = design_variable_update(design_variable,gradient,theta,kappa,options)

switch options.algorithm
    case 'level_set'
        
        phi_n = design_variable;   
        norm_g = sqrt(options.scalar_product(gradient,gradient));

        beta1 = sin((1-kappa)*theta)/sin(theta);
        beta2 = sin(kappa*theta)/sin(theta);
        phi = beta1*phi_n + beta2*gradient/norm_g;
               
        design_variable_new = phi;
        
    case 'Projected_gradient'
        gamma_n = design_variable;
        gamma_step = gamma_n - kappa*gradient;
        ub = options.ub(:);
        lb = options.lb(:);
        gamma = max(min(gamma_step,ub),lb);
        
        design_variable_new = gamma;
        
end