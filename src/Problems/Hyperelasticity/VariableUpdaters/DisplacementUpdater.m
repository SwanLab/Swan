classdef DisplacementUpdater < handle
    
    properties (Access = private)
        functional
        monitor

        tol
        maxIter

        pcgIterHistoryThisStep
        pcgMeanHistoryPerStep
    end

    methods (Access = public)

        function obj = DisplacementUpdater(cParams)
            obj.init(cParams);
        end

        function [u,rFun,costArray,iter] = update(obj,u,bc,costArray)

            obj.pcgIterHistoryThisStep = [];

            i = 0; err = 1; costOld = costArray(end);
            while (abs(err) > obj.tol) && (i < obj.maxIter)
                LHS = obj.functional.computeHessian(u);
                RHS = obj.functional.computeGradient(u,bc);
                u.setFValues(obj.computeDisplacement(LHS,RHS,u,bc));

                [err, cost] = obj.computeErrorCost(u,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;

                i = i+1;
                obj.monitor.printCost('iterU',i,cost,err);
                obj.monitor.update(length(costArray),{[],[cost],[],[]});
                obj.monitor.refresh(); 
            end
         
            if ~isempty(obj.pcgIterHistoryThisStep)
                meanCG_thisStep = mean(obj.pcgIterHistoryThisStep);
            else               
                meanCG_thisStep = NaN;
            end

            obj.pcgMeanHistoryPerStep(end+1) = meanCG_thisStep;

            figure(202); clf;
            bar(obj.pcgMeanHistoryPerStep); hold on; grid on;
            xlabel('Load step');
            ylabel('Mean PCG iterations per Newton');
            
            overallMean = mean(obj.pcgMeanHistoryPerStep,'omitnan');
            yline(overallMean, 'r', 'LineWidth', 2);

            title(sprintf('Mean PCG iters per Newton, per load step (overall mean = %.2f)', overallMean));
            legend('Mean per load step','Overall mean','Location','best');  

            rFun = obj.computeReactions(LHS,u,bc);
            iter = i;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.functional = cParams.functional;
            obj.monitor    = cParams.monitor;
            obj.tol        = cParams.tolerance;
            obj.maxIter    = cParams.maxIter;
            obj.pcgIterHistoryThisStep = [];
            obj.pcgMeanHistoryPerStep  = [];
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

        function xNew = updateWithNewton(obj,LHS,RHS,x)
    
    
            tol=10^(-8);
            maxIter = 1000000;
            Milu = obj.createILUpreconditioner(LHS);
            
            
            % LHSf = @(x) LHS*x;
            % 
            % [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x);
            
            [deltaX,flag,relres,iter,resvec] = pcg(LHS,-RHS,tol,maxIter,[],[]); 
            deltaX2 = -(LHS\RHS);

            obj.pcgIterHistoryThisStep(end+1) = iter;

            resSol = norm(deltaX - deltaX2) / norm(deltaX2); 
            resPcg = resvec / norm(-RHS);


            figure(101); clf;
            semilogy(resPcg,'LineWidth',1.5); hold on;
            semilogy(length(resPcg), resSol, 'o', 'MarkerSize',8, 'LineWidth',1.5);
            grid on;
            xlabel('iteration');
            ylabel('residual / error');
            legend('PCG relative residual ||b-Ax||/||b||','||deltaX - deltaX2|| / ||deltaX2||','Location','best');
            title(sprintf('pcg: flag=%d, relres=%.3e, iter=%d',flag,relres,iter));


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

        function Milu = createILUpreconditioner(~,LHS)

            s.LHS = sparse(LHS);
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);

        end

    end

end