function [compliance,gradient_compliance] = compliance_function(design_variable,epsilon_Le_kernel,kernel_case,equilibrium,diff_react_equation,ritz_gradient_representation,problembsc,Msmooth,h_C_0)

    global post_info
    [~,gamma_nodal,gamma_reg_gp] = diff_react_equation(design_variable,epsilon_Le_kernel,kernel_case,'gamma');
    [compliance,gradient_compliance,structural_values] = equilibrium(gamma_reg_gp,problembsc,h_C_0);
    
    gradient_compliance = diff_react_equation(gradient_compliance(:),epsilon_Le_kernel,kernel_case,'gradient');
    gradient_compliance = ritz_gradient_representation(gradient_compliance);
    gradient_compliance = Msmooth*gradient_compliance;
    
    post_info.gamma_reg_gp = gamma_reg_gp;
    post_info.structural_values = structural_values;
    post_info.gamma_nodal = gamma_nodal;
    post_info.compliance = compliance;
    post_info.gradient = gradient_compliance;
    post_info.xold = design_variable;
    post_info.Ch_star = problembsc.Ch_star;

end