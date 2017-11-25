classdef ShFunc_Compliance< Shape_Functional
    properties
        h_C_0=1; %compliance incial
    end
    methods
        function computef(obj,x,physicalProblem,interpolation,filter)  
            mass=filter.Msmooth;
            rho=filter.getP0fromP1(x);
            matProps=interpolation.computeMatProp(rho);
            %compute compliance
            compliance=physicalProblem.variables.d_u'*physicalProblem.RHS;  
            compliance=compliance/abs(obj.h_C_0);
            %compute gradient
            strain = physicalProblem.variables.strain;  
            gradient_compliance = zeros(physicalProblem.mesh.nelem,physicalProblem.geometry.ngaus); 
            for igaus=1:physicalProblem.geometry.ngaus
                for istre=1:physicalProblem.dim.nstre
                    for jstre = 1:physicalProblem.dim.nstre
                    gradient_compliance(:,igaus) = gradient_compliance(:,igaus) + (squeeze(-strain(igaus,istre,:)).*squeeze(matProps.dC(istre,jstre,:,igaus)).*squeeze(strain(igaus,jstre,:)));
                    end
                end
            end 
            gradient_compliance=gradient_compliance/abs(obj.h_C_0);
            gradient_compliance=filter.getP1fromP0(gradient_compliance);
            gradient_compliance = mass*gradient_compliance;
            
            obj.value=compliance;
            obj.gradient=gradient_compliance;
        end
    end
end
