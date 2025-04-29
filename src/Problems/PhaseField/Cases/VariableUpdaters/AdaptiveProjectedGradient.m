classdef AdaptiveProjectedGradient < handle

    properties (Access = private)
        functional
        monitor
        projectedGradient
        totIterPhi
    end
    
    methods (Access = public)
        function obj = AdaptiveProjectedGradient(cParams)
            obj.init(cParams);
        end

        function [phi,varargout] = update(obj,LHS,RHS,phi,varargin)
            u = varargin{1};
            bc = varargin{2};
            costOld = varargin{3};
            PG = obj.projectedGradient;

            phiNew = PG.update(RHS,phi);
            [err,~] = computeErrorCost(obj,u,phiNew,bc,costOld);
            while(err>0 && ~PG.isTooSmall())
                PG.decreaseStepLength();
                phiNew = PG.update(RHS,phi);
                [err,~] = computeErrorCost(obj,u,phiNew,bc,costOld);
            end
            phi = phiNew;
            tau = PG.tau;
            varargout{1} = tau;
            PG.increaseStepLength(10);
        end

        function updateBounds(obj,ub,lb)
            obj.projectedGradient.updateBounds(ub,lb)
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional        = cParams.functional;
            obj.monitor           = cParams.monitor;
            obj.projectedGradient = setProjectedGradient(obj,cParams);
            obj.totIterPhi        = 1;
        end

        function PG = setProjectedGradient(obj,cParams)
                s.ub = 1;
                s.lb = cParams.initPhi.fValues;
                s.tauMax = 1e10;
                s.tau = cParams.solver.tau;
                PG = ProjectedGradient(s);
        end

        function [e, cost] = computeErrorCost(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCost(u,phi,bc);
            e = cost - costOld;
        end

    end

end