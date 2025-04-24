classdef ProjectedNewton < handle

    properties (Access = private)
        upperBound
        lowerBound
    end
    
    methods (Access = public)

        function obj = ProjectedNewton(cParams)
            obj.init(cParams);
        end

        function var = update(obj,varargin)
            RHS = varargin{1}; 
            var = varargin{2}; 
            LHS = varargin{3};

            x  = var.fValues;
            xNew = obj.solve(LHS,RHS,x);
            xNew = obj.projectInBounds(xNew);
            var.setFValues(xNew); %% Change to designVariable
        end

        function updateBounds(obj,ub,lb)
            if isnumeric(ub)
                obj.upperBound = ub;
            else
                obj.upperBound = ub.fValues;
            end

            if isnumeric(lb)
                obj.lowerBound = lb;
            else
                obj.lowerBound = lb.fValues;
            end

        end

    end

    methods (Access = private)

        function init(obj,cParams)
            ub = 1;
            lb = cParams.initPhi;
            obj.updateBounds(ub,lb);
        end

        function xNew = solve(~,LHS,RHS,x)
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