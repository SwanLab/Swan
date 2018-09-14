classdef ShFunc_Compliance < Shape_Functional
    properties
        h_C_0; %compliance incial
        physProb
        interpolation
        rho
        matProps
    end
    
    methods
        function obj = ShFunc_Compliance(settings,postprocess_TopOpt)
            obj@Shape_Functional(settings);
            obj.physProb = FEM.create(settings.filename);
            obj.physProb.preProcess;
            obj.interpolation = Material_Interpolation.create(settings.TOL,settings.material,settings.method,settings.pdim);
            if settings.printing && settings.printing_physics
                obj.physProb.syncPostProcess(postprocess_TopOpt);
            end
        end
        function computeCostAndGradient(obj,x)
            obj.rho = obj.filter.getP0fromP1(x);
            obj.matProps = obj.interpolation.computeMatProp(obj.rho);
            obj.computeCostAndGradient_CORE;
        end
        function computeCostAndGradient_CORE(obj)
            % Compute compliance
            obj.physProb.setMatProps(obj.matProps);
            obj.physProb.computeVariables;
            
            compliance = obj.computeCompliance;
            
            % Compute gradient
            gradient_compliance = zeros(obj.physProb.geometry.interpolation.nelem,obj.physProb.element.quadrature.ngaus);
            for igaus = 1:obj.physProb.element.quadrature.ngaus
                for istre = 1:obj.physProb.element.nstre
                    for jstre = 1:obj.physProb.element.nstre
                        gradient_compliance(:,igaus) = gradient_compliance(:,igaus) + obj.updateGradient(igaus,istre,jstre);
                    end
                end
            end
            gradient_compliance = obj.filter.getP1fromP0(gradient_compliance);
            gradient_compliance = obj.Msmooth*gradient_compliance;
            
            if isempty(obj.h_C_0)
                obj.h_C_0 = compliance;
            end
            compliance = compliance/abs(obj.h_C_0);
            gradient_compliance = gradient_compliance/abs(obj.h_C_0);
            
            obj.value = compliance;
            obj.gradient = gradient_compliance;
        end
        
        function compliance = computeCompliance(obj)
            compliance = obj.physProb.variables.d_u'*obj.physProb.variables.fext;
        end
        
        function gradient_compliance = updateGradient(obj,igaus,istre,jstre)
            gradient_compliance = (squeeze(-obj.physProb.variables.strain(igaus,istre,:)).*squeeze(obj.matProps.dC(istre,jstre,:)).*squeeze(obj.physProb.variables.strain(igaus,jstre,:)));
        end
    end
end
