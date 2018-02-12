function [gamma_reg,lambda,cost_n,iter_level_set] = level_set_algorithm(problem)
global iter

lambda = problem.lambda0;
phi = problem.x0;
phi = phi/sqrt(problem.scalar_product(phi));
epsilon = problem.epsilon;
kernel_case = problem.kernel_case;

[cost,gradient,constraint,theta,norm_dif_rel,gamma_reg,gamma_reg_gp,structural_values,sub_cost_grad] = update_opt_variables_because_change_phi(phi,lambda,epsilon,problem,kernel_case);   
problem.print_function(iter,gamma_reg_gp,gamma_reg,gradient,phi,structural_values)
volume = problem.volume_from_constraint(constraint);
%plot_variables(problem.itert,0,cost,volume,lambda,theta,norm_dif_rel,0,1)
cost_n = cost; theta_n = theta; volume_n = volume; lambda_n = lambda; norm_dif_rel_n = norm_dif_rel; kappa_n = 0;gamma_old = gamma_reg; incre_gamma = 1;

stopping_criteria_opt = theta > problem.theta_min || constraint > problem.constraint_tol;
iter_level_set = 1;
while stopping_criteria_opt
    
    %Maximization thanks to lambda
    lambda = problem.update_lambda(lambda,constraint);
    [cost,gradient,theta] = update_opt_variables_because_change_lambda(phi,lambda,sub_cost_grad,epsilon,problem);
    
   % [cost_n,theta_n,volume_n,lambda_n] = update_saved_variables(cost_n,cost,theta_n,theta,volume_n,volume,lambda_n,lambda);   
    
    %Minimuzation thanks to phi (Line search)
    kappa = 1;stopping_Criteria_ls = 1; phi_ls = phi;k = 1; kfrac = 2;%1.04;
    while stopping_Criteria_ls
        phi_ls = slerp_update(phi,gradient,theta,kappa,kernel_case,problem);
        
        [cost_ls,gradient_ls,constraint_ls,theta_ls,norm_dif_rel_ls,gamma_reg_ls,gamma_reg_gp_ls,structural_values_ls,sub_cost_grad_ls] = update_opt_variables_because_change_phi(phi_ls,lambda,epsilon,problem,kernel_case);
        kappa = kappa/kfrac;
        incr_cost = (cost_ls - cost)/cost;
        stopping_Criteria_ls = incr_cost > 1e-6 && kappa >= problem.kappa_min; 
%         problem.print_function(k,gamma_reg_gp_ls,gamma_reg_ls,gradient_ls,phi_ls,structural_values_ls)
        cost_k(k) = cost_ls;
        vol_k(k) = problem.volume_from_constraint(constraint_ls);
        kappa_k(k) = kappa*kfrac;
        k = k + 1;

    end
    
    if kappa > problem.kappa_min
        phi = phi_ls; cost = cost_ls; gradient = gradient_ls;
        constraint = constraint_ls; theta = theta_ls; gamma_reg = gamma_reg_ls; norm_dif_rel = norm_dif_rel_ls;
        volume = problem.volume_from_constraint(constraint);
        gamma_reg_gp = gamma_reg_gp_ls; structural_values = structural_values_ls; sub_cost_grad = sub_cost_grad_ls;
    end
    
    
    incre_gamma(iter_level_set+1) = problem.scalar_product(gamma_reg - gamma_old)/problem.scalar_product(gamma_old);
    gamma_old = gamma_reg;
    [cost_n,theta_n,volume_n,lambda_n,norm_dif_rel_n,kappa_n] = update_saved_variables(cost_n,cost,theta_n,theta,volume_n,volume,lambda_n,lambda,norm_dif_rel_n,norm_dif_rel,kappa_n,kfrac*kappa);
    
    [stopping_criteria_opt,kernel_case] = compute_stopping_criteria(theta,constraint,problem,kernel_case);
    
    problem.print_function(iter,gamma_reg_gp,gamma_reg,gradient,phi,structural_values)

    plot_variables(problem.itert + iter_level_set,iter_level_set,cost_n,volume_n,lambda_n,theta_n,norm_dif_rel_n,kappa_n,incre_gamma)
    
    theta*180/pi
    volume
    
    minmax = [min(gamma_reg_gp),max(gamma_reg_gp)]
    
    kappa
    
    problem.scalar_product(gradient)
    
    iter_level_set = iter_level_set + 1;
    
end

end