classdef Optimizer < handle
    
    properties (Access = protected)
        designVariable
        dualVariable
        cost
        constraint
        outputFunction
        maxIter
        nIter
        targetParameters
        dualUpdater
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

        function createDualUpdater(obj,cParams)
            f               = DualUpdaterFactory();
            obj.dualUpdater = f.create(cParams);      
        end
        
        function isAcceptable = checkConstraint(obj)
            switch obj.constraintCase{1}
                case {'EQUALITY'}
                    isAcceptable = obj.checkEqualityConstraint();
                case {'INEQUALITY'}
                    isAcceptable = obj.checkInequalityConstraint();
            end
        end
        
    end

    methods (Access = private)

        function c = checkInequalityConstraint(obj)
            g = obj.constraint.value;
            c = g < obj.targetParameters.constr_tol;
        end

        function c = checkEqualityConstraint(obj)
            g = obj.constraint.value;
            c = abs(g) < obj.targetParameters.constr_tol;
        end

    end
    
end