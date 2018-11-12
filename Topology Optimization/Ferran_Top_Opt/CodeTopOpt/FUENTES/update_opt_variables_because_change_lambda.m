function [cost,gradient,theta] = update_opt_variables_because_change_lambda(design_variable,lambda_volume,penalty_volume,lambda_Perimeter,penalty_Perimeter,sub_cost_grad,epsilon,problem)

compliance = sub_cost_grad.compliance;
gradient_compliance = sub_cost_grad.gradient_compliance;
constraint_volume = sub_cost_grad.constraint_volume;
gradient_constraint_volume = sub_cost_grad.gradient_constraint_volume;
constraint_Perimeter = sub_cost_grad.constraint_Perimeter;
gradient_constraint_Perimeter = sub_cost_grad.gradient_constraint_Perimeter;


switch problem.volume_integration

    case 'natural_integration'
      [cost_mechanical,gradient_mechanical] = problem.assamble_cost_gradient_mechanical(compliance,gradient_compliance,constraint_volume,gradient_constraint_volume,constraint_Perimeter,gradient_constraint_Perimeter,...
                                                        lambda_volume,penalty_volume,lambda_Perimeter,penalty_Perimeter);
                                                    
      [cost,gradient] = problem.assamble_cost_gradient(cost_mechanical,gradient_mechanical,constraint_Perimeter,gradient_constraint_Perimeter,...
                                                        lambda_Perimeter,penalty_Perimeter);
      

    case 'regularized_integration'
        [cost_mechanical,gradient_mechanical] = problem.assamble_cost_gradient_mechanical(compliance,gradient_compliance,constraint_volume,gradient_constraint_volume,...
                                                        lambda_volume,penalty_volume);
        gradient_mechanical = problem.diff_react_equation(gradient_mechanical,epsilon,problem.kernel_case,'gradient');
        
        [cost,gradient] = problem.assamble_cost_gradient(cost_mechanical,gradient_mechanical,constraint_Perimeter,gradient_constraint_Perimeter,...
                                                        lambda_Perimeter,penalty_Perimeter);
      
        
        

end


gradient = problem.ritz_gradient_representation(gradient);                  


theta = problem.call_theta(design_variable,gradient); 
end