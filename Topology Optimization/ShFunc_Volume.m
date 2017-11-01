classdef ShFunc_Volume< Shape_Functional
    properties 
        Vfrac
    end
    methods 
        function obj=ShFunc_Volume(Vfrac)
            obj.Vfrac=Vfrac;
        end
        function [volume, gradient_volume]=computef(gamma,physicalProblem,interpolation)
            %Update phys problem
            matProps=interpolation.computeMatProp(gamma);
            physicalProblem.element.material.setMatProps(matProps)
            physicalProblem.computeVariables;
            M0 = sparse(1:dim.nelem,1:dim.nelem,dvolu);
            
            geometric_volume = sum(physicalProblem.element.computeMass(:));
            volume = sum(M0*gamma_gp);
            
            
            %volume = volume/geometric_volume;
            volume = volume/(geometric_volume*obj.Vfrac) - 1;
            
            gradient_volume = 1/(geometric_volume*obj.Vfrac);
            gradient_volume = gradient_volume*ones(size(element.conectivities,1),1);
            % constraint_gradient = M0*constraint_gradient;
            %           gradient_volume = diff_react_equation(gradient_volume,epsilon_Le_kernel,kernel_case,'gradient');
            %           gradient_volume = ritz_gradient_representation(gradient_volume);
            %           gradient_volume = Msmooth*gradient_volume;
        end
    end
end
