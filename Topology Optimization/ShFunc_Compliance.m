classdef ShFunc_Compliance < Shape_Functional
    properties
        h_C_0; %compliance incial
        physicalProblem
        interpolation
        rho
        matProps
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
            obj.interpolation = Material_Interpolation.create(settings.TOL,settings.material,settings.method);
        end
        
        function computef(obj,x)
            obj.rho = obj.filter.getP0fromP1(x);
            obj.matProps = obj.interpolation.computeMatProp(obj.rho);
            obj.computef_CORE;
        end
        
        function computef_CORE(obj)
            % Compute compliance
            obj.physicalProblem.setMatProps(obj.matProps);
            obj.physicalProblem.computeVariables;
            
            compliance = obj.computeCompliance;
            
            % Compute gradient
            gradient_compliance = zeros(obj.physicalProblem.mesh.nelem,obj.physicalProblem.geometry.ngaus);
            for igaus = 1:obj.physicalProblem.geometry.ngaus
                for istre = 1:obj.physicalProblem.dim.nstre
                    for jstre = 1:obj.physicalProblem.dim.nstre
                        gradient_compliance(:,igaus) = gradient_compliance(:,igaus) + obj.updateGradient(igaus,istre,jstre);
                    end
                end
            end
            
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
        
        function compliance = computeCompliance(obj)
            compliance = obj.physicalProblem.variables.d_u'*obj.physicalProblem.variables.fext;
        end
        
        function gradient_compliance = updateGradient(obj,igaus,istre,jstre)
            gradient_compliance = (squeeze(-obj.physicalProblem.variables.strain(igaus,istre,:)).*squeeze(obj.matProps.dC(istre,jstre,:)).*squeeze(obj.physicalProblem.variables.strain(igaus,jstre,:)));
        end
    end
end
