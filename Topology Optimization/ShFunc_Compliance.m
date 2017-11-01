classdef ShFunc_Compliance< Shape_Functional
    properties 
    end
    methods (Static)
        function [compliance, gradient_compliance]=computef(gamma,physicalProblem,interpolation)
            %Update phys problem
            matProps=interpolation.computeMatProp(gamma);
            physicalProblem.element.material.setMatProps(matProps);
            physicalProblem.computeVariables;
            
            compliance=physicalProblem.physicalProperties.U'*physicalProblem.physicalProperties.fext;    
            gradient_compliance = P_operator'*M0*gradient_compliance;
            %commented:L2, uncommented:h1
            %gradient_compliance = (epsilon_P0_gradient^2*Stiff_smooth + Msmooth)\(Msmooth*gradient);
            gradient_compliance = physicalProblem.element.computeMass*gradient_compliance;
    
        end
    end
end
