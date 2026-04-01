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

        function [phi,varargout] = update(obj,hessian,gradient,phi,varargin)
            u = varargin{1};
            bc = varargin{2};
            costOld = varargin{3};
            PG = obj.projectedGradient;

            phi.updateOld();
            phi = PG.update(gradient,phi,bc.phi);
            [err,~] = computeErrorCost(obj,u,phi,bc.u,costOld);
            while(err>0 && ~PG.isTooSmall())
                phi.recoverOld();
                PG.decreaseStepLength();
                phi = PG.update(gradient,phi,bc.phi);
                [err,~] = computeErrorCost(obj,u,phi,bc.u,costOld);
            end
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
                s.lb = cParams.initPhi.fun.fValues;
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