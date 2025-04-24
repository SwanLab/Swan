classdef PhaseFieldDamageUpdater < handle
    
    properties (Access = private)
        functional
        monitor
        tol
        maxIter
        solver 
    end

    methods (Access = public)

        function obj = PhaseFieldDamageUpdater(cParams)
            obj.init(cParams);
        end

        function [phi,costArray,iter] = update(obj,u,phi,bc,costArray)
            iter = 0; err = 1; costOld = costArray(end);
            while (abs(err) > obj.tol) && (iter < obj.maxIter)
                LHS = obj.functional.computePhaseFieldLHS(u,phi);
                RHS = obj.functional.computePhaseFieldRHS(u,phi);
                phi = obj.solver.update(RHS,phi,LHS);

                [err, cost] = obj.computeErrorCost(u,phi,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;

                iter = iter+1;
                obj.monitor.printCost('iterPhi',iter,cost,err);
                obj.monitor.update(length(costArray),{[],[],[],[],[cost],[],[]});
                obj.monitor.refresh();
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
            switch cParams.solverType
                case 'Gradient'
                    obj.solver = ProjectedGradient(cParams);
                case 'Newton'
                    obj.solver = ProjectedNewton(cParams);
            end
        end

        function [e, cost] = computeErrorCost(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCostFunctional(u,phi,bc);
            e = cost - costOld;
        end

    end

end