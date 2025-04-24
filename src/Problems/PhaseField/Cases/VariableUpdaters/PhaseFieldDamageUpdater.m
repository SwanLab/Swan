classdef PhaseFieldDamageUpdater < OptimizerPhaseField
    
    properties (Access = private)
        functional
        tol
        maxIter
        solver 
    end

    methods (Access = public)

        function obj = PhaseFieldDamageUpdater(cParams)
            obj.init(cParams);
        end

        function [phi,costArray,iter] = update(u,phi,bc,costArray)
            iter = 0; err = 0; costOld = costArray(end);
            while (abs(err) > obj.tol) && (iter < obj.maxIter)
                LHS = obj.functional.computePhaseFieldLHS(u,phi);
                RHS = obj.functional.computePhaseFieldRHS(u,phi);
                phi = obj.solver.update(RHS,phi,LHS);

                [err, cost] = computeErrorCostFunctional(u,phi,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;

                iter = iter+1;
                obj.monitor.printCost('iterPhi',iter,cost,err);
                obj.monitor.update(length(costArray),{[],[],[],[],[cost],[],[]});
                obj.monitor.refresh();
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.tol        = cParams.toleranceDamage;
            obj.maxIter    = cParams.maxIterDamage;
            switch cParams.solverType
                case 'Gradient'
                    obj.solver = ProjectedGradient(cParams);
                case 'Newton'
                    obj.solver = ProjectedNewton(cParams);
            end
        end

        function updateBounds(obj,ub,lb)
            obj.solver.updateBounds(ub,lb);
        end

        function xNew = updateWithGradient(~,RHS,x,tau)
            deltaX = -tau.*RHS;
            xNew = x + deltaX; 
        end

        function [e, cost] = computeErrorCost(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCostFunctional(u,phi,bc);
            e = cost - costOld;
        end

        function printCost(obj,name,iter,cost,e)
            if obj.print
                X = sprintf('%s:%d / cost: %.8e  (diff:%.8e) \n',name,iter,cost,e);
                fprintf(X);
            end
        end

    end

end