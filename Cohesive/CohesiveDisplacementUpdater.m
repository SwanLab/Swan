classdef CohesiveDisplacementUpdater < handle
    
    properties (Access = private)
        functional
        monitor

        tol
        maxIter
    end

    methods (Access = public)

        function obj = CohesiveDisplacementUpdater(cParams)
            obj.init(cParams);
        end

        function [u,rFun,resNormArray,iter] = update(obj,u,bc,resNormArray)
            i = 0;
            resNorm = inf;
        
            while (resNorm > obj.tol) && (i < obj.maxIter)
        
                [LHSSec,LHSTan] = obj.functional.computeHessian(u);
                RHS = obj.functional.computeGradient(u,bc);
                u.setFValues(obj.computeDisplacement(LHSTan,RHS,u,bc));
        
        
                RHS = obj.functional.computeGradient(u,bc);
        
                [~,RHSRed] = fullToReduced(obj,LHSTan,RHS,bc);
                resNorm = norm(RHSRed);
        
                resNormArray(end+1) = resNorm;
        
                i = i + 1;
        
                fprintf('Iteration %d: Reduced Residual Norm = %.4e\n', ...
                    i, resNorm);
            end
        
            rFun = obj.computeReactions(LHSSec,u,bc);
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

        function rFun = computeReactions(~,LHS,u,bc)
            constrainedDofs = bc.dirichlet_dofs;
            uVec = reshape(u.fValues',[u.nDofs 1]);
            KR   = LHS(constrainedDofs,:);
            rVec = zeros(size(uVec));
            rVec(constrainedDofs,1) = -KR*uVec;
            
            s.fValues = reshape(rVec,[flip(size(u.fValues))])';
            s.order   = 'P1';
            s.mesh    = u.mesh;
            rFun = LagrangianFunction(s);
        end

    end

end