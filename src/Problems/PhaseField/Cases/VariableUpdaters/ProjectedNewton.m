classdef ProjectedNewton < handle

    properties (Access = private)
        upperBound
        lowerBound
    end
    
    methods (Access = public)

        function obj = ProjectedNewton(cParams)
            obj.init(cParams);
        end

        function var = update(~,varargin)
            RHS = varargin{1}; 
            var = varargin{2}; 
            LHS = varargin{3};

            x  = var.fValues;
            xNew = solve(LHS,RHS,x);
            xNew  = projectInBounds(xNew);
            var.update(xNew);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            ub = cParams.ub;
            lb = cParams.lb;
            obj.updateBounds(ub,lb);
        end

        function updateBounds(obj,ub,lb)
            obj.upperBound = ub;
            obj.lowerBound = lb;
        end 

        function xNew = solve(LHS,RHS,x)
            deltaX = -LHS\RHS;
            xNew = x + deltaX;
        end

        function xP = projectInBounds(obj,x)
            xLB = obj.lowerBound;
            xUB = obj.upperBound;
            xP = min(max(xLB, x),xUB);
        end

    end

end