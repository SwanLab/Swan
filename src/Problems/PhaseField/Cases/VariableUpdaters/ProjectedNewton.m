classdef ProjectedNewton < handle

    properties (Access = private)
        upperBound
        lowerBound
    end
    
    methods (Access = public)

        function obj = ProjectedNewton(cParams)
            obj.init(cParams);
        end

        function [phi,varargout] = update(obj,hessian,gradient,phi,varargin)
            x  = phi.fun.fValues;
            xNew = obj.solve(hessian,gradient,x);
            xNew = obj.projectInBounds(xNew);
            phi.update(xNew);
            varargout{1} = [];
        end
        
        function updateBounds(obj,ub,lb)
            if ~isnumeric(ub)
                obj.upperBound = ub.fValues;
            else
                obj.upperBound = ub;
            end
            if ~isnumeric(lb)
                obj.lowerBound = lb.fValues;
            else
                obj.lowerBound = lb;
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            ub = 1;
            lb = cParams.initPhi.fun;
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