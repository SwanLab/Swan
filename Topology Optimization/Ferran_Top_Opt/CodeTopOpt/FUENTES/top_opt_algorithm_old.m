function [gamma_reg,lambda_volume,cost_n,iter_level_set] = top_opt_algorithm(problem)
global iter

lambda_volume = problem.lambda_volume_0;
penalty_volume = problem.penalty_volume_0;
lambda_Perimeter = problem.lambda_Perimeter_0;
penalty_Perimeter = problem.penalty_Perimeter_0;
design_variable = problem.x0;
epsilon_Le_kernel = problem.epsilon_Le_kernel;
epsilon_iter_perimeter = problem.epsilon_iter_perimeter;
kernel_case = problem.kernel_case;

%plot_variables(problem.itert,0,cost,volume,lambda,theta,norm_dif_rel,0,1)
cost_n = []; theta_n = []; volume_n = []; lambda_volume_n = []; Perimeter_n = []; kappa_n = 0; gamma_old = []; incre_gamma_n = [];

%stopping_criteria_opt = theta > problem.theta_min || constraint_volume > problem.constraint_tol;
gamma_old = rand(size(problem.gamma0));
incre_gamma = 1;
iter_level_set = 1;
problem.itert = 0;

for iepsilon = 1:length(epsilon_iter_perimeter)
    epsilon_Perimeter = epsilon_iter_perimeter(iepsilon);
    
    [cost,gradient,constraint_volume,constraint_Perimeter,theta,norm_dif_rel,gamma_reg,gamma_reg_gp,structural_values,sub_cost_grad] = update_opt_variables_because_change_design_variable(design_variable,lambda_volume,penalty_volume,lambda_Perimeter,penalty_Perimeter,epsilon_Le_kernel,epsilon_Perimeter,problem,kernel_case);
    volume = problem.volume_from_constraint(constraint_volume);
    Perimeter = problem.Perimeter_from_constraint(constraint_Perimeter);
    %incre_gamma = 1;%sqrt(problem.scalar_product(gamma_reg - gamma_old,gamma_reg - gamma_old))/sqrt(problem.scalar_product(gamma_old,gamma_old));
    gamma_old = gamma_reg;
    [cost_n,theta_n,volume_n,lambda_volume_n,Perimeter_n,kappa_n,incre_gamma_n] = update_saved_variables(cost_n,cost,theta_n,theta,volume_n,volume,lambda_volume_n,lambda_volume,Perimeter_n,Perimeter,kappa_n,0,incre_gamma_n,incre_gamma);
   
    [stopping_criteria_opt,kernel_case] = compute_stopping_criteria(theta,incre_gamma,constraint_volume,problem,kernel_case,problem.algorithm,iepsilon, problem.penalty_volume_0);
    problem.print_function(iter,gamma_reg_gp,gamma_reg,gradient,design_variable,structural_values)
    plot_variables(problem.itert + iter_level_set,iter_level_set,cost_n,volume_n,lambda_volume_n,theta_n,Perimeter_n,kappa_n,incre_gamma_n)
        
    
    %stopping_criteria_opt = 1;
    
    while stopping_criteria_opt
        
        %Maximization thanks to lambda
        penalty_volume = problem.update_penalty_volume(penalty_volume);
        lambda_volume = problem.update_lambda(lambda_volume,penalty_volume,constraint_volume);
        penalty_Perimeter = problem.update_penalty_Perimeter(penalty_Perimeter);
        lambda_Perimeter = problem.update_lambda(lambda_Perimeter,penalty_Perimeter,constraint_Perimeter);
       
        [cost,gradient,theta] = update_opt_variables_because_change_lambda(design_variable,lambda_volume,penalty_volume,lambda_Perimeter,penalty_Perimeter,sub_cost_grad,epsilon_Le_kernel,problem);
        
         
        % [cost_n,theta_n,volume_n,lambda_n] = update_saved_variables(cost_n,cost,theta_n,theta,volume_n,volume,lambda_n,lambda);
        
        %Minimuzation thanks to design_variable (Line search)
        kappa = inicialize_kappa(design_variable,gradient,problem);
        stopping_Criteria_ls = 1; design_variable_ls = design_variable;k = 1; kfrac = 2;%1.04;
        while stopping_Criteria_ls
            design_variable_ls = design_variable_update(design_variable,gradient,theta,kappa,kernel_case,problem);
            
            [cost_ls,gradient_ls,constraint_volume_ls,constraint_Perimeter_ls,theta_ls,norm_dif_rel_ls,gamma_reg_ls,gamma_reg_gp_ls,structural_values_ls,sub_cost_grad_ls] = update_opt_variables_because_change_design_variable(design_variable_ls,lambda_volume,penalty_volume,lambda_Perimeter,penalty_Perimeter,epsilon_Le_kernel,epsilon_Perimeter,problem,kernel_case);
            kappa = kappa/kfrac;
            incr_cost = (cost_ls - cost)/abs(cost);
            stopping_Criteria_ls = incr_cost > 0 && kappa >= problem.kappa_min;
            %         problem.print_function(k,gamma_reg_gp_ls,gamma_reg_ls,gradient_ls,design_variable_ls,structural_values_ls)
            cost_k(k) = cost_ls;
            vol_k(k) = problem.volume_from_constraint(constraint_volume_ls);
            kappa_k(k) = kappa*kfrac;
            k = k + 1;
            
        end
        
        if kappa > problem.kappa_min
            design_variable = design_variable_ls; cost = cost_ls; gradient = gradient_ls;
            constraint_volume = constraint_volume_ls; constraint_Perimeter = constraint_Perimeter_ls;
            theta = theta_ls; gamma_reg = gamma_reg_ls; norm_dif_rel = norm_dif_rel_ls;
            volume = problem.volume_from_constraint(constraint_volume); Perimeter = problem.Perimeter_from_constraint(constraint_Perimeter);
            gamma_reg_gp = gamma_reg_gp_ls; structural_values = structural_values_ls; sub_cost_grad = sub_cost_grad_ls;
        end
        
        
        incre_gamma = sqrt(problem.scalar_product(gamma_reg - gamma_old,gamma_reg - gamma_old))/sqrt(problem.scalar_product(gamma_old,gamma_old))
        gamma_old = gamma_reg;
        [cost_n,theta_n,volume_n,lambda_volume_n,Perimeter_n,kappa_n,incre_gamma_n] = update_saved_variables(cost_n,cost,theta_n,theta,volume_n,volume,lambda_volume_n,lambda_volume,Perimeter_n,Perimeter,kappa_n,kfrac*kappa,incre_gamma_n,incre_gamma);
        
        [stopping_criteria_opt,kernel_case] = compute_stopping_criteria(theta,incre_gamma,constraint_volume,problem,kernel_case,problem.algorithm,iepsilon, problem.penalty_volume_0);
        
        problem.print_function(iter,gamma_reg_gp,gamma_reg,gradient,design_variable,structural_values)
        
        plot_variables(problem.itert + iter_level_set,iter_level_set,cost_n,volume_n,lambda_volume_n,theta_n,Perimeter_n,kappa_n,incre_gamma_n)
        
        theta*180/pi
        volume
        minmax = [min(gamma_reg_gp),max(gamma_reg_gp)]
        kappa
        sqrt(problem.scalar_product(gamma_reg,gamma_reg))
        iepsilon
        
        iter_level_set = iter_level_set + 1;
        
    end
    
end

end