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
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        type
    end
    
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f    = OptimizerFactory();
            obj  = f.create(cParams);
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
            obj.outputFunction    = cParams.outputFunction.monitoring;
        end
        
    end
    
end