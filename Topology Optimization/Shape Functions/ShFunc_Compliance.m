classdef ShFunc_Compliance < ShFunWithElasticPdes
    
    methods (Access = public)
        
        function obj = ShFunc_Compliance(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            obj.init(cParams);
            fileName = cParams.femSettings.fileName;
            obj.createEquilibriumProblem(fileName);
        end
        
    end
    
    methods (Access = protected)
        
        function solvePDEs(obj)
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C)
            obj.physicalProblem.computeVariables();
        end
        
        function computeFunctionValue(obj)
            u = obj.physicalProblem.variables.d_u;
            f = obj.physicalProblem.variables.fext;
            obj.value = f'*u;
        end
        
        function g = updateGradient(obj,igaus,istre,jstre)
            e    = obj.physicalProblem.variables.strain;
            ei   = squeeze(e(igaus,istre,:));
            ej   = squeeze(e(igaus,jstre,:));
            g = zeros(length(ei),1,obj.nVariables);
            for ivar = 1:obj.nVariables
                dCij = squeeze(obj.homogenizedVariablesComputer.dC(istre,jstre,ivar,:));
                g(:,1,ivar) = -ei.*dCij.*ej;
            end
        end
        
    end
    
end

