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
            bc = varargin{2};
            fDofs = bc.phi.free_dofs;
            [LHS,RHS] = obj.fullToReduced(hessian,gradient,bc.phi);

            xFree  = phi.fun.fValues(fDofs);
            xNew   = phi.fun.fValues;
            xNew(fDofs) = obj.solve(LHS,RHS,xFree);
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

        function [LHS,RHS] = fullToReduced(~,LHS,RHS,bc)
            free_dofs = bc.free_dofs;
            LHS = LHS(free_dofs, free_dofs);
            RHS = RHS(free_dofs);
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