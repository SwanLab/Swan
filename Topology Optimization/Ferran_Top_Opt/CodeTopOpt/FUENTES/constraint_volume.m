function [constraint,constraint_gradient] = constraint_volume(design_variable,Vfrac,dim,element,problembsc,coordinates,Msmooth,M0,algorithm,volume_integration,epsilon_Le_kernel,kernel_case,diff_react_equation,ritz_gradient_representation)

global post_info
[~,~,gamma_gp] = diff_react_equation(design_variable,epsilon_Le_kernel,kernel_case,'gamma');

switch volume_integration
    case 'natural_integration'
        switch algorithm
            case 'level_set'
                phi = design_variable;
                [~,~,~,~,volume,geometric_volume] = cal_vol_mat(phi,dim,element,problembsc,coordinates);
                
            case 'Projected_gradient'
                gamma = design_variable;
                geometric_volume = sum(Msmooth(:));
                volume = sum(Msmooth)*gamma;
                
        end
        
    case 'regularized_integration'
          geometric_volume = sum(Msmooth(:));
          volume = sum(M0*gamma_gp);
        
end
post_info.volume = volume/geometric_volume;
constraint = volume/(geometric_volume*Vfrac) - 1;

constraint_gradient = 1/(geometric_volume*Vfrac);
constraint_gradient = constraint_gradient*ones(size(element.conectivities,1),1);
% constraint_gradient = M0*constraint_gradient;
constraint_gradient = diff_react_equation(constraint_gradient,epsilon_Le_kernel,kernel_case,'gradient');
constraint_gradient = ritz_gradient_representation(constraint_gradient);
constraint_gradient = Msmooth*constraint_gradient;
end