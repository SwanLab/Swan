classdef ShFunc_NonSelfAdjoint_Compliance < Shape_Functional
    properties
        h_C_0; %compliance incial
        forces_adjoint
    end
    methods
        function obj=ShFunc_NonSelfAdjoint_Compliance(settings)
            obj.forces_adjoint=Preprocess.getBC_adjoint(settings.filename);
%            obj@Shape_Functional(settings);
        end
        function computef(obj,x,physicalProblem,interpolation,filter)  
            rho=filter.getP0fromP1(x);
            matProps=interpolation.computeMatProp(rho);
            physicalProblem.setMatProps(matProps);
            physicalProblem.computeVariables;
            
            RHS=physicalProblem.RHS;
            strain = physicalProblem.variables.strain;
            
            adjointProblem=physicalProblem;
            adjointProblem.bc.neunodes=obj.forces_adjoint;
            adjointProblem.bc.neunodes(:,3)=-adjointProblem.bc.neunodes(:,3);
            adjointProblem.preProcess;
            adjointProblem.computeVariables;
            strain_adjoint=adjointProblem.variables.strain;
            
            compliance=adjointProblem.variables.d_u'*RHS;
            
            %compute gradient            
            gradient_compliance = zeros(physicalProblem.mesh.nelem,physicalProblem.geometry.ngaus); 
            for igaus=1:physicalProblem.geometry.ngaus
                for istre=1:physicalProblem.dim.nstre
                    for jstre = 1:physicalProblem.dim.nstre
                    gradient_compliance(:,igaus) = gradient_compliance(:,igaus) + (squeeze(strain(igaus,istre,:)).*squeeze(matProps.dC(istre,jstre,:)).*squeeze(strain_adjoint(igaus,jstre,:)));
                    end
                end
            end 
            
            gradient_compliance=filter.getP1fromP0(gradient_compliance);
            gradient_compliance = filter.Msmooth*gradient_compliance;
            if isempty(obj.h_C_0)
                obj.h_C_0=compliance;
            else
                compliance=compliance/abs(obj.h_C_0);
                gradient_compliance=gradient_compliance/abs(obj.h_C_0);
            end            
            obj.value=compliance;
            obj.gradient=gradient_compliance;
        end
    end
end
