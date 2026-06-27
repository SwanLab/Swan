classdef DisplacementUpdater < handle
    
    properties (Access = private)
        functional
        monitor

        tol
        maxIter
    end

    methods (Access = public)

        function obj = DisplacementUpdater(cParams)
            obj.init(cParams);
        end

        function [u,F,costArray,iter] = update(obj,u,bc,costArray)
            i = 0; err = 1; costOld = costArray(end); 
            % normOld = inf; useSecant = 0;

            while (abs(err) > obj.tol) && (i < obj.maxIter)
                [LHSSec,LHSTan] = obj.functional.computeHessian(u);
                RHS = obj.functional.computeGradient(u,bc);

                [~,RHSRed] = fullToReduced(obj,LHSTan,RHS,bc);
                normRHS = norm(RHSRed);


                u.setFValues(obj.computeDisplacement(LHSSec,RHS,u,bc));
                
                [err, cost] = obj.computeErrorCost(u,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;
                
                i = i+1;

                normOld = normRHS;
                fprintf('Iteration %d: Reduced Residual Norm = %.4e\n', i, normRHS);
            end

            obj.functional.updateDamageOld(u);
            F = obj.computeForceVector(LHSSec,u);
            iter = i;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional = cParams.functional;
            obj.monitor    = cParams.monitor;
            obj.tol        = cParams.tolerance;
            obj.maxIter    = cParams.maxIter;
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

        function [e, cost] = computeErrorCost(obj,u,bc,costOld)
            cost = obj.functional.computeCost(u,bc); % To include extWork
            e = (cost - costOld)/cost;
        end

        function F = computeForceVector(~,LHS,u)
            uVec = reshape(u.fValues',[u.nDofs 1]);
            F = LHS*uVec;
        end

    end

end