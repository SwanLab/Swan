classdef Optimizer < handle

    properties (Access = public)
        nIter = 0
        outputFunction
    end
    
    properties (Access = protected)
        designVariable
        dualVariable
        cost
        constraint
        %outputFunction
        maxIter 
        targetParameters
        dualUpdater
        primalUpdater
        constraintCase
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        type
    end
    
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f   = OptimizerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function initOptimizer(obj,cParams)
            obj.nIter             = 0;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.maxIter           = cParams.maxIter;
            obj.targetParameters  = cParams.targetParameters;
            obj.constraintCase    = cParams.constraintCase;
            obj.outputFunction    = cParams.outputFunction.monitoring;
        end

        function createPrimalUpdater(obj,cParams)
            f                 = PrimalUpdaterFactory();
            obj.primalUpdater = f.create(cParams);
        end

        function createDualUpdater(obj,cParams)
            f               = DualUpdaterFactory();
            obj.dualUpdater = f.create(cParams);      
        end
        
        function isAcceptable = checkConstraint(obj)
            for i = 1:length(obj.constraint.value)
                switch obj.constraintCase{i}
                    case {'EQUALITY'}
                        fine = obj.checkEqualityConstraint(i);
                    case {'INEQUALITY'}
                        fine = obj.checkInequalityConstraint(i);
                end
                if fine
                    isAcceptable = true;
                else
                    isAcceptable = false;
                end
            end
        end
        
    end

    methods (Access = private)

        function c = checkInequalityConstraint(obj,i)
            g = obj.constraint.value(i);
            c = g < obj.targetParameters.constr_tol;
        end

        function c = checkEqualityConstraint(obj,i)
            g = obj.constraint.value(i);
            c = abs(g) < obj.targetParameters.constr_tol;
        end

    end
    
end