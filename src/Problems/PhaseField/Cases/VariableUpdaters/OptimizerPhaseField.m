classdef OptimizerPhaseField < handle
    
    properties (Access = private)
        functional
        monitor

        tol
        maxIter
    end

    properties (Access = private)
        displacementUpdater
        damageUpdater
    end

    methods (Access = public)

        function obj = OptimizerPhaseField(cParams)
            obj.init(cParams);
        end

        function [u,phi,F,costArray,iterMax] = compute(obj,u,phi,bc,costArray)
            % Staggered scheme here
            iterMax.u = 1; iterMax.phi = 1; iterMax.stag = 1;

            iter = 0; err = 1; costOld = costArray(end);
            while (abs(err) > obj.tol) && (iter < obj.maxIter)
                [u,F,costArray,iterU]   = obj.updateDisplacement(u,phi,bc,costArray);
                iterMax.u = obj.checkIterMax(iterU,iterMax.u);

                [phi,costArray,iterPhi] = obj.updateDamage(obj,u,phi,bc,costArray);
                iterMax.phi = obj.checkIterMax(iterPhi,iterMax.phi);

                [err, cost] = obj.computeErrorCost(u,phi,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;
        
                iter = iter+1;
                obj.monitor.printCost('iterStag',iter,cost,err);
                obj.monitor.update(length(costArray),{[],[],[],[],[cost],[],[]});
                obj.monitor.refresh();
            end
        end

    end

    methods (Access = protected)

        function [e, cost] = computeErrorCost(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCostFunctional(u,phi,bc);
            e = cost - costOld;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional = cParams.functional;
            obj.monitor    = cParams.monitor;
            obj.displacementUpdater = PhaseFieldDisplacementUpdater(cParams);
            obj.damageUpdater       = PhaseFieldDamageUpdater(cParams);
        end

        function [u,F,costArray,iter] = updateDisplacement(obj,u,phi,bc,costArray)
            dispUpdater = obj.displacementUpdater;
            [u,F,costArray,iter] = dispUpdater.update(obj,u,phi,bc,costArray);
        end

        function [phi,costArray,iter] = updateDamage(obj,u,phi,bc,costArray)
            dmgUpdater = obj.damageUpdater;
            [phi,costArray,iter] = dmgUpdater.update(obj,u,phi,bc,costArray);
        end

        function iterMax = checkIterMax(obj,iter,iterMax)
            if iter > iterMax
                iterMax = iter;
            end
        end

    end

end