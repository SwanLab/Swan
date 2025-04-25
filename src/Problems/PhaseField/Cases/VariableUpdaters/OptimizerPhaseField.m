classdef OptimizerPhaseField < handle
    
    properties (Access = protected)
        functional
        monitor
    end

    properties (Access = private)
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

        function [u,phi,F,costArray,iter] = compute(obj,u,phi,bc,costArray)
            iter.u = 1; iter.phi = 1; iter.stag = 1;
            i = 0; err = 1; costOld = costArray(end);
            while (abs(err) > obj.tol) && (i < obj.maxIter)
                [u,F,costArray,iterU]   = obj.updateDisplacement(u,phi,bc,costArray);
                iter.u = max(iterU,iter.u);

                [phi,costArray,iterPhi] = obj.updateDamage(u,phi,bc,costArray);
                iter.phi = max(iterPhi,iter.phi);

                [err, cost] = obj.computeErrorCost(u,phi,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;
        
                i = i+1;
                obj.monitor.printCost('iterStag',i,cost,err);
                obj.monitor.update(length(costArray),{[],[],[],[],[cost],[],[]});
                obj.monitor.refresh();
            end
            iter.stag = i;
            obj.damageUpdater.updateBounds(1,phi);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional = cParams.functional;
            obj.monitor    = cParams.monitor;
            obj.tol        = cParams.tolerance.stag;
            obj.maxIter    = cParams.maxIter.stag;
            obj.displacementUpdater = PhaseFieldDisplacementUpdater(cParams);
            obj.damageUpdater       = PhaseFieldDamageUpdater(cParams);
        end

        function [u,F,costArray,iter] = updateDisplacement(obj,u,phi,bc,costArray)
            dispUpdater = obj.displacementUpdater;
            [u,F,costArray,iter] = dispUpdater.update(u,phi,bc,costArray);
        end

        function [phi,costArray,iter] = updateDamage(obj,u,phi,bc,costArray)
            dmgUpdater = obj.damageUpdater;
            [phi,costArray,iter] = dmgUpdater.update(u,phi,bc,costArray);
        end

        function [e, cost] = computeErrorCost(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCost(u,phi,bc);
            e = cost - costOld;
        end

    end

end