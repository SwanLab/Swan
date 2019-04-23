classdef ShFunc_Compliance < ShFunWithElasticPdes
   
    methods (Access = public)
        
        function obj = ShFunc_Compliance(cParams)
            obj.init(cParams);
            obj.createEquilibriumProblem(cParams.filename);
            obj.createHomogenizedVariablesComputer(cParams);            
        end
        
    end
    
    methods (Access = protected)
        
        function solvePDEs(obj)
            obj.physProb.setC(obj.homogenizedVariablesComputer.C)
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
            dCij = squeeze(obj.homogenizedVariablesComputer.dC(istre,jstre,:));
            g    = -ei.*dCij.*ej;
        end
        
    end
    
end

