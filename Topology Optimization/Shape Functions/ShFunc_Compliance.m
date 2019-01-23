classdef ShFunc_Compliance < ShFunWithElasticPdes
   
    methods (Access = public)
        
        function obj = ShFunc_Compliance(settings)
            obj@ShFunWithElasticPdes(settings);
            obj.createEquilibriumProblem(settings.filename);
        end
        
    end
    
    methods (Access = protected)
        
        function solvePDEs(obj)
            obj.physProb.setMatProps(obj.matProps);
            obj.physProb.computeVariables();
        end
        
        function computeFunctionValue(obj)
            u = obj.physProb.variables.d_u;
            f = obj.physProb.variables.fext;
            obj.value = f'*u;
        end
        
        function g = updateGradient(obj,igaus,istre,jstre)
            e    = obj.physProb.variables.strain;
            ei   = squeeze(e(igaus,istre,:));
            ej   = squeeze(e(igaus,jstre,:));
            dCij = squeeze(obj.matProps.dC(istre,jstre,:));
            g    = -ei.*dCij.*ej;
        end
        
    end
    
end

