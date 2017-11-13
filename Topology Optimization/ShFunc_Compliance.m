classdef ShFunc_Compliance< Shape_Functional
    properties 
    end
    methods (Static)
        function [compliance, gradient_compliance]=computef(rho,physicalProblem,interpolation)
            %Update phys problem
            matProps=interpolation.computeMatProp(rho);
            physicalProblem.setMatProps(matProps);
            physicalProblem.computeVariables;
            
            compliance=physicalProblem.variables.d_u'*physicalProblem.RHS;    
            
            strain = physicalProblem.variables.strain;
            stress = physicalProblem.variables.stress;
            
            %derivada C
            P_interp = zeros(physicalProblem.dim.nstre,physicalProblem.dim.nstre,length(rho));
            P_interp(1,1,:) = matProps.c1 + matProps.c2;
            P_interp(2,2,:) = matProps.c1 + matProps.c2;
            P_interp(1,2,:) = matProps.c2;
            P_interp(2,1,:) = matProps.c2;
            P_interp(3,3,:) = matProps.c1;            
            
            gradient_compliance = zeros(physicalProblem.mesh.nelem,physicalProblem.geometry.ngaus);
            
            for igaus=1:physicalProblem.geometry.ngaus
                for istre=1:physicalProblem.dim.nstre
                    for jstre = 1:physicalProblem.dim.nstre
                        gradient_compliance(:,igaus) = gradient_compliance(:,igaus) + (squeeze(stress(istre,:,igaus))*squeeze(P_interp(istre,jstre,:))*squeeze(strain(jstre,:,igaus)))';
                    end
                end
            end

            
            %gradient_compliance = P_operator'*M0*gradient_compliance;
            %commented:L2, uncommented:h1
            %gradient_compliance = (epsilon_P0_gradient^2*Stiff_smooth + Msmooth)\(Msmooth*gradient);
            %gradient_compliance =physicalProblem.element.computeMass*gradient_compliance;
    
        end
    end
end
