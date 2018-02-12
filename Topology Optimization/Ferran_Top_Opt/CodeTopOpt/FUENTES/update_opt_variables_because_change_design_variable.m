function [cost,gradient,constraint_volume,constraint_Perimeter,theta,norm_dif_rel,gamma_nodal,gamma_reg_gp,structural_values,sub_cost_grad] = update_opt_variables_because_change_design_variable(design_variable,lambda_volume,penalty_volume,lambda_Perimeter,penalty_Perimeter,epsilon_Le_kernel,epsilon_Perimeter,problem,kernel_case)

[~,gamma_nodal,gamma_reg_gp] = problem.diff_react_equation(design_variable,epsilon_Le_kernel,kernel_case,'gamma');
%gamma_reg_gp = problem.nodal_2_elem(gamma_reg);
[compliance,gradient_compliance,structural_values] = problem.objective(gamma_reg_gp);

switch problem.volume_integration

    case 'natural_integration'
        [constraint_volume,gradient_constraint_volume,constr_eq_volume] = problem.constraint_volume(design_variable,gamma_reg_gp,lambda_volume,penalty_volume);
        [constraint_Perimeter,gradient_constraint_Perimeter,constr_eq_perimeter] = problem.constraint_Perimeter(design_variable,gamma_reg_gp,epsilon_Perimeter,lambda_Perimeter,penalty_Perimeter);
        gradient_compliance = problem.diff_react_equation(gradient_compliance,epsilon_Le_kernel,kernel_case,'gradient');
        [cost_mechanical,gradient_mechanical] = problem.assamble_cost_gradient_mechanical(compliance,gradient_compliance,constraint_volume,gradient_constraint_volume,...
                                                        lambda_volume,penalty_volume);
        
        [cost,gradient] = problem.assamble_cost_gradient(cost_mechanical,gradient_mechanical,constraint_Perimeter,gradient_constraint_Perimeter,...
                                                        lambda_Perimeter,penalty_Perimeter);
      
                                                    
        sub_cost_grad.compliance = compliance;
        sub_cost_grad.gradient_compliance = gradient_compliance;
        sub_cost_grad.constraint_volume = constraint_volume;
        sub_cost_grad.gradient_constraint_volume = gradient_constraint_volume;
        sub_cost_grad.constraint_Perimeter = constraint_Perimeter;
        sub_cost_grad.gradient_constraint_volume = gradient_constraint_volume;
        sub_cost_grad.constr_eq_volume = constr_eq_volume;
        sub_cost_grad.constr_eq_perimeter = constr_eq_perimeter;

    case 'regularized_integration'
        [constraint_volume,gradient_constraint_volume,constr_eq_volume] = problem.constraint_volume(design_variable,gamma_reg_gp,lambda_volume,penalty_volume);
        [constraint_Perimeter,gradient_constraint_Perimeter,constr_eq_perimeter] = problem.constraint_Perimeter(design_variable,gamma_reg_gp,epsilon_Perimeter,lambda_Perimeter,penalty_Perimeter);
        [cost_mechanical,gradient_mechanical] = problem.assamble_cost_gradient_mechanical(compliance,gradient_compliance,constraint_volume,gradient_constraint_volume,...
                                                            lambda_volume,penalty_volume);
        gradient_mechanical = problem.diff_react_equation(gradient_mechanical,epsilon_Le_kernel,kernel_case,'gradient');
        
        [cost,gradient] = problem.assamble_cost_gradient(cost_mechanical,gradient_mechanical,constraint_Perimeter,gradient_constraint_Perimeter,...
                                                        lambda_Perimeter,penalty_Perimeter);
        
        sub_cost_grad.compliance = compliance;
        sub_cost_grad.gradient_compliance = gradient_compliance;
        sub_cost_grad.constraint_volume = constraint_volume;
        sub_cost_grad.gradient_constraint_volume = gradient_constraint_volume;
        sub_cost_grad.constraint_Perimeter = constraint_Perimeter;
        sub_cost_grad.gradient_constraint_Perimeter = gradient_constraint_Perimeter;
        sub_cost_grad.constr_eq_volume = constr_eq_volume;
        sub_cost_grad.constr_eq_perimeter = constr_eq_perimeter;

end

gradient = problem.ritz_gradient_representation(gradient);


[theta,norm_dif_rel] = problem.call_theta(design_variable,gradient); 

end