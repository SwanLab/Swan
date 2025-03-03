classdef DualUpdater_AugmentedLagrangian < handle
    
    properties (Access = private)
       dualVariable
       penalty
       constraint
       nConstr
       constraintCase
    end
        
    methods (Access = public)
        
        function obj = DualUpdater_AugmentedLagrangian(cParams)
            obj.init(cParams)
        end
        
        function update(obj)
            obj.updateDual();
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
            obj.nConstr        = length(cParams.constraintCase);
        end

        function compute(obj,i)
            l   = obj.dualVariable.fun.fValues(i);
            rho = obj.penalty;
            c   = obj.constraint.value(i);
            l   = l + rho.*c;
            obj.dualVariable.fun.fValues(i) = l;
        end

        function updateDual(obj)
            for i = 1:obj.nConstr
                switch obj.constraintCase{i}
                    case 'INEQUALITY'
                        isZero = obj.checkDual(i);
                        if isZero
                            obj.dualVariable.fun.fValues(i) = 0;
                        else
                            obj.compute(i);
                            obj.dualVariable.fun.fValues(i) = max(0,obj.dualVariable.fun.fValues(i));
                        end
                    otherwise
                        obj.compute(i);
                end
            end
        end

        function isZero = checkDual(obj,i)
            g      = obj.constraint.value(i,1);
            l      = obj.dualVariable.fun.fValues(i,1);
            rho    = obj.penalty;
            isZero = g < -l/rho;
        end

    end
    
  
end