function [constraint,constraint_gradient] = constraint_Perimeter(design_variable,epsilon,algorithm,Msmooth,M0,Perimeter_target,element,coordinates,diff_react_equation,ritz_gradient_representation)

global post_info
gamma_reg = diff_react_equation(design_variable,epsilon,'PDE','gamma');
switch algorithm
    case 'level_set'
        rhs = faireF2(coordinates',element.conectivities',design_variable);
      
    otherwise % case 'Projected_gradient'
        rhs = Msmooth*design_variable;
        
end

Perimeter = 0.5/epsilon*((1 - gamma_reg)'*rhs);
Perimeter_gradient = 0.5/epsilon*(1 - 2*gamma_reg);

constraint = Perimeter/Perimeter_target - 1;
constraint_gradient = Perimeter_gradient/Perimeter_target;
constraint_gradient = ritz_gradient_representation(constraint_gradient);
constraint_gradient = Msmooth*constraint_gradient;

post_info.perimeter = Perimeter;

end