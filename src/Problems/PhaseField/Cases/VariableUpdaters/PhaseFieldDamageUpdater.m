classdef PhaseFieldDamageUpdater < OptimizerPhaseField
    
    properties (Access = private)
        functional
        tol
        maxIter
        solver 

        monitor
        print
    end

    methods (Access = public)

        function obj = PhaseFieldDamageUpdater(cParams)
            obj.init(cParams);
        end

        function [u,F,costArray,iter] = update(u,phi,bc,costArray)
            iter = 0; err = 0; costOld = costArray(end);
            while (abs(err) > obj.tol) && (iter < obj.maxIter)
                LHS = obj.functional.computePhaseFieldLHS(u,phi);
                RHS = obj.functional.computePhaseFieldRHS(u,phi);
                phi = obj.solver.update(RHS,phi,LHS);

                [err, cost] = obj.computeErrorCostFunctional(u,phi,bc,costOld);
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
            obj.functional = cParams.functional;
            obj.tol        = cParams.tolerance;
            obj.maxIter    = cParams.maxIter;
            obj.monitor    = cParams.monitor;
            obj.print      = cParams.print;
            switch cParams.solverType
                case 'Gradient'
                    obj.solver = ProjectedGradient(cParams);
                case 'Newton'
                    obj.solver = ProjectedNewton(cParams);
            end
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