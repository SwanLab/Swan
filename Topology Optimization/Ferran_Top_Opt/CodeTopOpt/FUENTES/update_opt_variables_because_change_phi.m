function [cost,gradient,constraint,theta,norm_dif_rel,gamma_nodal,gamma_reg_gp,structural_values,sub_cost_grad] = update_opt_variables_because_change_design_variable(design_variable,lambda,epsilon,problem,kernel_case)

[~,gamma_nodal,gamma_reg_gp] = problem.diff_react_equation(design_variable,epsilon,kernel_case,'gamma',design_variable);
%gamma_reg_gp = problem.nodal_2_elem(gamma_reg);
[compliance,gradient_compliance_p,structural_values] = problem.objective(gamma_reg_gp);
gradient_compliance = problem.diff_react_equation(gradient_compliance_p,epsilon,kernel_case,'gradient',design_variable);
%gradient_compliance2 = problem.diff_react_equation((problem.A_dvolu)*gradient_compliance_p,epsilon,'gradient_boundary',design_variable);
[constraint,gradient_constraint] = problem.constraint(design_variable);
[cost,gradient] = problem.assamble_cost_gradient(compliance,gradient_compliance,constraint,gradient_constraint,lambda);
%gradient = problem.diff_react_equation((problem.A_nodal_2_gauss)'*gradient_p,epsilon,'P1');

[theta,norm_dif_rel] = problem.call_theta(design_variable,gradient); 

sub_cost_grad.compliance = compliance;
sub_cost_grad.gradient_compliance = gradient_compliance;
sub_cost_grad.constraint = constraint;
sub_cost_grad.gradient_constraint = gradient_constraint;

end