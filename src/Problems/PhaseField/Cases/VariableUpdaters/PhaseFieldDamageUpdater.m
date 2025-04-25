classdef PhaseFieldDamageUpdater < handle
    
    properties (Access = private)
        functional
        monitor
        tol
        maxIter
        solver
        totIter
    end

    methods (Access = public)

        function obj = PhaseFieldDamageUpdater(cParams)
            obj.init(cParams);
        end

        function [phi,costArray,iter] = update(obj,u,phi,bc,costArray)
            iter = 1; err = 1; costOld = costArray(end);
            while (abs(err) > obj.tol) && (iter < obj.maxIter)
                LHS = obj.functional.computePhaseFieldLHS(u,phi);
                RHS = obj.functional.computePhaseFieldRHS(u,phi);
                [phi,tau] = obj.solver.update(LHS,RHS,phi,u,bc,costOld);

                [err, cost] = obj.computeErrorCost(u,phi,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;

                obj.monitor.printCost('iterPhi',iter,cost,err);
                obj.monitor.update(obj.totIter,{[],[],[],[],[],[],[tau]})
                obj.monitor.updateAndRefresh(length(costArray),{[],[],[],[],[cost],[],[]});
                iter = iter+1;
                obj.totIter = obj.totIter + 1;
            end
        end

        function updateBounds(obj,ub,lb)
            obj.solver.updateBounds(ub,lb);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional = cParams.functional;
            obj.monitor    = cParams.monitor;
            obj.tol        = cParams.tolerance.phi;
            obj.maxIter    = cParams.maxIter.phi;
            obj.totIter    = 1;
            switch cParams.solverType
                case 'Gradient'
                    obj.solver = AdaptiveProjectedGradient(cParams);
                case 'Newton'
                    obj.solver = ProjectedNewton(cParams);
            end
        end

        function [e, cost] = computeErrorCost(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCost(u,phi,bc);
            e = cost - costOld;
        end

    end

end