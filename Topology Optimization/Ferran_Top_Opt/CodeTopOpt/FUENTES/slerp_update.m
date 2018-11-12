function phi = design_variable_update(phi_n,gradient,theta,kappa,kernel_case,problem)

switch problem.Algorithm
    case 'level_set'

norm_g = sqrt(problem.scalar_product(gradient));
norm_phi =sqrt(problem.scalar_product(phi_n));
switch kernel_case
    
    case 'P1_kernel'
        %%slerp
        
        beta1 = sin((1-kappa)*theta)/sin(theta);
        beta2 = sin(kappa*theta)/sin(theta);
        
        phi = beta1*phi_n + beta2*gradient/norm_g;
        
        %%gradient
    case 'P0_kernel'
       % phi = phi_n + kappa*norm_phi/norm_g*gradient;
       % phi = phi/sqrt(problem.scalar_product(phi));
        beta1 = sin((1-kappa)*theta)/sin(theta);
        beta2 = sin(kappa*theta)/sin(theta);
        
        phi = beta1*phi_n + beta2*gradient/norm_g;
       
        
    case 'Le_epsilon'
        beta1 = sin((1-kappa)*theta)/sin(theta);
        beta2 = sin(kappa*theta)/sin(theta);
        
        phi = beta1*phi_n + beta2*gradient/norm_g;
        
        
       
end

    case 'Projected_Gradient'
        
        

end