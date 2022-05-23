classdef DualUpdater_AugmentedLagrangian < handle
    
    properties (Access = private)
       dualVariable
       penalty
       constraint
       nConstr
       constraintCase
       index
    end
        
    methods (Access = public)
        
        function obj = DualUpdater_AugmentedLagrangian(cParams)
            obj.init(cParams)
        end
        
        function update(obj)
            switch obj.constraintCase{1}
                case {'EQUALITY'}
                    obj.updateDual();
                case {'INEQUALITY'}
                    obj.updateInequalityDual();
            end
        end

        function updatePenalty(obj,rho)
            obj.penalty = rho;
        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.constraint     = cParams.constraint;
            obj.constraintCase = cParams.constraintCase;
            obj.dualVariable   = cParams.dualVariable;
            obj.nConstr        = cParams.constraint.nSF;
        end

        function updateDual(obj)
            l   = obj.dualVariable.value;
            rho = obj.penalty;
            c   = obj.constraint.value;
            l   = l + rho.*c;
            obj.dualVariable.value = l;
        end

        function updateInequalityDual(obj)
            if obj.isNotZero()
                obj.updateDual();
            else
                obj.dualVariable.value = zeros(obj.index,1);
            end
        end

        function c = isNotZero(obj)
            g     = obj.constraint.value;
            l     = obj.dualVariable.value;
            rho   = obj.penalty;
            c     = g + l/rho > 0;
            index = find(c == 0) 
        end

    end
    
  
end