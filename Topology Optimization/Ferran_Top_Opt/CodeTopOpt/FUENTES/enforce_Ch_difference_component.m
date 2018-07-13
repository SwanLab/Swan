function [Ch_diff2,gradient_Ch_diff] = enforce_Ch_difference_component ( design_variable,compliancefun,epsilon_Le_kernel,kernel_case,diff_react_equation,ritz_gradient_representation,Msmooth,typeCdifference,component,initial_values )

switch typeCdifference
    case 'isotropy'
        [costval,gradientval] = isotropy_equation2D (design_variable,compliancefun);
        %[costval,gradientval] = isotropy_equation2D_scaled (design_variable,compliancefun);
    case 'CCstar'
        [costval,gradientval] = CCstar_equation (design_variable,compliancefun);
    otherwise
        error('Type C difference not detected.');
end
Ch_diff2 = (costval(component));
gradient_Ch_diff = gradientval(:,component);

Ch_diff2 = Ch_diff2/initial_values(component);
gradient_Ch_diff = gradient_Ch_diff/initial_values(component);

gradient_Ch_diff = diff_react_equation(gradient_Ch_diff,epsilon_Le_kernel,kernel_case,'gradient');
gradient_Ch_diff = ritz_gradient_representation(gradient_Ch_diff);
gradient_Ch_diff = Msmooth*gradient_Ch_diff;

end

