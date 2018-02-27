classdef ShFunc_Compliance < Shape_Functional
    properties
        h_C_0; %compliance incial
        physicalProblem
    end
    methods
        function obj = ShFunc_Compliance(settings)
            obj@Shape_Functional(settings);
            switch settings.ptype
                case 'MACRO'
                    obj.physicalProblem = Physical_Problem(settings.filename);
                case 'MICRO'
                    obj.physicalProblem = Physical_Problem_Micro(settings.filename);
            end
            obj.physicalProblem.preProcess;
        end
        
%         function obj = preProcess(obj)
%             rho = obj.filter.getP0fromP1(obj.x);
%             matProps = obj.interpolation.computeMatProp(rho);
%             obj.physicalProblem.setMatProps(matProps);
%         end
        
        function computef(obj,x,interpolation)
            rho = obj.filter.getP0fromP1(x);
            matProps = interpolation.computeMatProp(rho);
            
            %compute compliance
            obj.physicalProblem.setMatProps(matProps);
            obj.physicalProblem.computeVariables;
            
            %compliance = obj.physicalProblem.variables.d_u'*obj.physicalProblem.RHS;
            compliance = obj.physicalProblem.variables.d_u'*obj.physicalProblem.variables.fext;
            
            %compute gradient
            strain = obj.physicalProblem.variables.strain;
            stress = obj.physicalProblem.variables.stress;
            fobj = 0;
            gradient_compliance = zeros(obj.physicalProblem.mesh.nelem,obj.physicalProblem.geometry.ngaus);
            for igaus = 1:obj.physicalProblem.geometry.ngaus
                for istre = 1:obj.physicalProblem.dim.nstre
                    for jstre = 1:obj.physicalProblem.dim.nstre
                        gradient_compliance(:,igaus) = gradient_compliance(:,igaus) + (squeeze(-strain(igaus,istre,:)).*squeeze(matProps.dC(istre,jstre,:)).*squeeze(strain(igaus,jstre,:)));
                    end
                    fobj = fobj + (squeeze(strain(igaus,istre,:)).*squeeze(stress(igaus,istre,:)))'*obj.physicalProblem.geometry.dvolu(:,igaus);
                end
            end
            
            %% !! NOTE: INVERSE ORDER THAN NON-SELF-ADJOINT COMPLIANCE, MISTAKE? !!
            
            gradient_compliance = obj.filter.getP1fromP0(gradient_compliance);
            gradient_compliance = obj.filter.Msmooth*gradient_compliance;
            
            if isempty(obj.h_C_0)
                obj.h_C_0 = compliance;
            else
                compliance = compliance/abs(obj.h_C_0);
                gradient_compliance = gradient_compliance/abs(obj.h_C_0);
            end
            
            obj.value = compliance;
            obj.gradient = gradient_compliance;
        end
    end
end
