classdef Optimizer_PG < Optimizer_Unconstrained
    
    properties  (GetAccess = public, SetAccess = private)
        optimality_tol
    end
    
    properties (GetAccess = public, SetAccess = protected)
        type = 'PROJECTED GRADIENT'
    end
    
    properties (Access = private)
       upperBound
       lowerBound
    end
    
    methods (Access = public)
        
        function obj = Optimizer_PG(cParams)
            obj@Optimizer_Unconstrained(cParams);
            obj.upperBound = cParams.ub;
            obj.lowerBound = cParams.lb;
        end
        
        function compute(obj)
            x_n      = obj.designVariable.value;
            gradient = obj.objectiveFunction.gradient;
            x_new = x_n-obj.lineSearch.value*gradient;
            ub = obj.upperBound*ones(length(x_n(:,1)),1);
            x2 = reshape(x_n,[],2); norm(x2(:,1)-x2(:,2))
            lb = obj.lowerBound*ones(length(x_n(:,1)),1);
            x_new = max(min(x_new,ub),lb);
            obj.designVariable.update(x_new);
            l2Norm = obj.designVariable.computeL2normIncrement();
            obj.optimalityCond = sqrt(l2Norm);
        end
        
    end
    
end