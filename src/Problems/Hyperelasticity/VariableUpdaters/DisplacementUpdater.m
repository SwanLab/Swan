classdef DisplacementUpdater < handle

    properties (Access = private)
        functional
        monitor

        tol
        maxIter
        useLineSearch
        residualTol
        viscosity
    end

    methods (Access = public)

        function obj = DisplacementUpdater(cParams)
            obj.init(cParams);
        end

        function [u,rFun,costArray,iter] = update(obj,u,bc,costArray)
            if obj.useLineSearch
                [u,LHS,costArray,iter] = obj.updateWithProximalNewton(u,bc,costArray);
            else
                [u,LHS,costArray,iter] = obj.updateWithFullNewton(u,bc,costArray);
            end
            rFun = obj.computeReactions(LHS,u,bc);
        end

    end

    methods (Static, Access = public)

        function [x,tau] = solveWithPositiveShift(K,b)
            n   = size(K,1);
            x   = nan(n,1);
            tau = 0;
            if DisplacementUpdater.isFinite(K,b)
                [R,q,tau] = DisplacementUpdater.factorizeWithPositiveShift(K);
                x(q) = R\(R'\b(q));
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional = cParams.functional;
            obj.monitor    = cParams.monitor;
            obj.tol        = cParams.tolerance;
            obj.maxIter    = cParams.maxIter;
            obj.useLineSearch = isfield(cParams,'lineSearch') && cParams.lineSearch;
            if obj.useLineSearch
                obj.residualTol = cParams.residualTolerance;
                obj.viscosity   = cParams.viscosity;
            end
        end

        function [u,LHS,costArray,i] = updateWithFullNewton(obj,u,bc,costArray)
            i = 0; err = 1; costOld = costArray(end);
            while (abs(err) > obj.tol) && (i < obj.maxIter)
                LHS = obj.functional.computeHessian(u);
                RHS = obj.functional.computeGradient(u,bc);
                u.setFValues(obj.computeDisplacement(LHS,RHS,u,bc));
                [err, cost] = obj.computeErrorCost(u,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;
                i = i+1;
                obj.printIteration(i,cost,err,costArray);
            end
        end

        function [u,LHS,costArray,i] = updateWithProximalNewton(obj,u,bc,costArray)
            free = bc.free_dofs;
            uVec = obj.toDofVector(u);
            eta  = [];
            i = 0;
            while i < obj.maxIter
                LHS = obj.functional.computeHessian(u);
                RHS = obj.functional.computeGradient(u,bc);
                g   = RHS(free);
                if norm(g) < obj.residualTol
                    break
                end
                K = LHS(free,free);
                D = obj.computeProximalWeights(K);
                if isempty(eta)
                    eta = obj.viscosity;
                end
                [uVec,cost,alpha,eta] = obj.takeDampedStep(K,g,uVec,free,u,bc,eta,D);
                if isempty(alpha)
                    i = obj.maxIter;
                    break
                end
                u.setFValues(obj.toNodalValues(uVec,u));
                i = i+1;
                obj.printIteration(i,cost,norm(g),costArray);
            end
            costArray(end+1) = obj.computeTrialCost(uVec,u,bc);
        end

        function [uVec,cost,alpha,eta] = takeDampedStep(obj,K,g,uVec,free,u,bc,eta,D)
            uRef = uVec(free);
            n    = numel(free);
            while eta < 1e6
                costFun = @(v) obj.computeTrialCost(v,u,bc) + eta/2*sum(D.*(v(free)-uRef).^2);
                [uTrial,cost,alpha] = obj.takeLineSearchStep(K + eta*spdiags(D,0,n,n),g,uVec,free,costFun);
                if ~isempty(alpha)
                    uVec = uTrial;
                    if alpha == 1
                        eta = eta/10;
                    else
                        eta = 10*eta;
                    end
                    return
                end
                eta = max(100*eta,obj.viscosity);
            end
            alpha = [];
            cost  = [];
        end

        function [uVec,cost,alpha] = takeLineSearchStep(obj,K,g,uVec,free,costFun)
            cost0 = costFun(uVec);
            [dNewton,tau] = DisplacementUpdater.solveWithPositiveShift(K,g);
            dNewton = -dNewton;
            [xNew,cost,alpha] = obj.searchAlong(dNewton,g'*dNewton,cost0,uVec,free,costFun);
            if tau > 0
                d = obj.computeNegativeCurvatureDirection(K,g,norm(dNewton),tau);
                [xCurv,costCurv,alphaCurv] = obj.searchAlong(d,g'*d,cost0,uVec,free,costFun);
                if ~isempty(alphaCurv) && (isempty(alpha) || costCurv < cost)
                    xNew  = xCurv;
                    cost  = costCurv;
                    alpha = alphaCurv;
                end
            end
            uVec(free) = xNew;
        end

        function [xNew,cost,acceptedAlpha] = searchAlong(~,d,slope,cost0,uVec,free,costFun)
            alpha = 1;
            acceptedAlpha = [];
            xNew  = uVec(free);
            cost  = cost0;
            for k = 1:30
                uTrial       = uVec;
                uTrial(free) = uVec(free) + alpha*d;
                trialCost = costFun(uTrial);
                if isreal(trialCost) && isfinite(trialCost) && trialCost <= cost0 + 1e-4*alpha*min(slope,0) + 1e-12*abs(cost0)
                    xNew  = uTrial(free);
                    cost  = trialCost;
                    acceptedAlpha = alpha;
                    break
                end
                alpha = alpha/2;
            end
        end

        function d = computeNegativeCurvatureDirection(~,K,g,len,tau)
            [v,lambda] = eigs((K+K')/2,1,-2*tau);
            d = zeros(size(g));
            if lambda < 0
                s = -sign(g'*v);
                if s == 0
                    s = 1;
                end
                d = s*len*v/norm(v);
            end
        end

        function printIteration(obj,i,cost,err,costArray)
            obj.monitor.printCost('iterU',i,cost,err);
            obj.monitor.update(length(costArray),{[],[cost],[],[]});
            obj.monitor.refresh();
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

        function cost = computeTrialCost(obj,uVec,uIn,bc)
            trial = copy(uIn);
            trial.setFValues(obj.toNodalValues(uVec,uIn));
            cost = obj.functional.computeCost(trial,bc);
        end

        function v = toDofVector(~,fun)
            v = reshape(fun.fValues',[fun.nDofs 1]);
        end

        function f = toNodalValues(~,v,fun)
            f = reshape(v,[flip(size(fun.fValues))])';
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
            s.order   = u.order;
            s.mesh    = u.mesh;
            rFun = LagrangianFunction(s);
        end

    end

    methods (Static, Access = private)

        function D = computeProximalWeights(K)
            D = abs(full(diag(K)));
            D = max(D,1e-12*mean(D));
        end

        function [R,q,tau] = factorizeWithPositiveShift(K)
            n       = size(K,1);
            minimal = max(1e-8*mean(abs(diag(K))),eps);
            tau     = 0;
            [R,p,q] = chol(K,'vector');
            while p > 0
                tau     = max(10*tau,minimal);
                [R,p,q] = chol(K + tau*speye(n),'vector');
            end
            if tau > 0
                [R,~,q] = chol(K + 4*tau*speye(n),'vector');
            end
        end

        function itIs = isFinite(K,b)
            itIs = all(isfinite(nonzeros(K))) && all(isfinite(b));
        end

    end

end
