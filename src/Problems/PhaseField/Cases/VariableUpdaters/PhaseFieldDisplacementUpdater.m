classdef PhaseFieldDisplacementUpdater < handle
    
    properties (Access = private)
        functional
        monitor

        tol
        maxIter
    end

    methods (Access = public)

        function obj = PhaseFieldDisplacementUpdater(cParams)
            obj.init(cParams);
        end

        function [u,F,costArray,iter] = update(obj,u,phi,bc,costArray)
            i = 0; err = 1; costOld = costArray(end);
            while (abs(err) > obj.tol) && (i < obj.maxIter)
                LHS = obj.functional.computeElasticLHS(u,phi);
                RHS = obj.functional.computeElasticRHS(u,phi,bc);
                u.setFValues(obj.computeDisplacement(LHS,RHS,u,bc));

                [err, cost] = obj.computeErrorCost(u,phi,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;

                i = i+1;
                obj.monitor.printCost('iterU',i,cost,err);
                obj.monitor.update(length(costArray),{[],[],[],[],[cost],[],[]});
                obj.monitor.refresh();
                
            end
            F = obj.computeForceVector(LHS,u);
            iter = i;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional = cParams.functional;
            obj.monitor    = cParams.monitor;
            obj.tol        = cParams.tolerance.u;
            obj.maxIter    = cParams.maxIter.u;
        end

        function uOut = computeDisplacement(obj,LHSfull, RHSfull,uIn,bc)
            [LHS,RHS] = fullToReduced(obj,LHSfull,RHSfull,bc);
            if ~isempty(LHS)
                uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
                uOutVec = uInVec;

                uInFree = uInVec(bc.free_dofs);
                uOutFree = obj.updateWithNewton(LHS,RHS,uInFree);
                uOutVec(bc.free_dofs) = uOutFree;
                uOut = reshape(uOutVec,[flip(size(uIn.fValues))])';
            else
                uOut = uIn.fValues;
            end
        end

        function [LHS,RHS] = fullToReduced(~,LHS,RHS,bc)
            free_dofs = bc.free_dofs;
            LHS = LHS(free_dofs, free_dofs);
            RHS = RHS(free_dofs);
        end

        function xNew = updateWithNewton(~,LHS,RHS,x)
            deltaX = -LHS\RHS;
            xNew = x + deltaX;
        end

        function [e, cost] = computeErrorCost(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCostFunctional(u,phi,bc);
            e = cost - costOld;
        end

        function F = computeForceVector(~,LHS,u)
            uVec = reshape(u.fValues',[u.nDofs 1]);
            F = LHS*uVec;
        end

    end

end