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
        angleUpdater
        damageUpdater
    end

    methods (Access = public)

        function obj = OptimizerPhaseField(cParams)
            obj.init(cParams);
        end

        function [u,theta,phi,F,costArray,iter] = compute(obj,u,theta,phi,bc,costArray)
            iter.u = 1; iter.phi = 1; iter.stag = 1;
            i = 0; err = 1; costOld = costArray(end);
            while (abs(err) > obj.tol) && (i < obj.maxIter)
                [u,F,costArray,iterU]   = obj.updateDisplacement(u,theta,phi,bc.u,costArray);
                iter.u = max(iterU,iter.u);

                [theta] = obj.updateOrientation(u,theta,phi);

                [phi,costArray,iterPhi] = obj.updateDamage(u,theta,phi,bc,costArray);
                iter.phi = max(iterPhi,iter.phi);

                [err, cost] = obj.computeErrorCost(u,theta,phi,bc.u,costOld);
                costArray(end+1) = cost;
                costOld = cost;
        
                i = i+1;
                obj.monitor.printCost('iterStag',i,cost,err);
                obj.monitor.update(length(costArray),{[],[],[],[],[cost],[],[]});
                obj.monitor.refresh();
            end
            iter.stag = i;
            obj.damageUpdater.updateBounds(1,phi.fun);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional = cParams.functional;
            obj.monitor    = cParams.monitor;
            obj.tol        = cParams.tolerance.stag;
            obj.maxIter    = cParams.maxIter.stag;
            obj.displacementUpdater = PhaseFieldDisplacementUpdater(cParams);
            obj.angleUpdater        = PhaseFieldAngleUpdater();
            obj.damageUpdater       = PhaseFieldDamageUpdater(cParams);
        end

        function [u,F,costArray,iter] = updateDisplacement(obj,u,theta,phi,bc,costArray)
            dispUpdater = obj.displacementUpdater;
            [u,F,costArray,iter] = dispUpdater.update(u,theta,phi,bc,costArray);
        end

        function [theta] = updateOrientation(obj,u,theta,phi)
            thetaUpdater = obj.angleUpdater;
            [theta] = thetaUpdater.update(u,theta,phi);
        end

        function [phi,costArray,iter] = updateDamage(obj,u,theta,phi,bc,costArray)
            dmgUpdater = obj.damageUpdater;
            [phi,costArray,iter] = dmgUpdater.update(u,theta,phi,bc,costArray);
        end

        function [e, cost] = computeErrorCost(obj,u,theta,phi,bc,costOld)
            cost = obj.functional.computeCost(u,theta,phi,bc);
            e = cost - costOld;
        end

    end

end