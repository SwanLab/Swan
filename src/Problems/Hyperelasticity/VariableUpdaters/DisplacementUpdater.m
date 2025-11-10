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

        function [u,rFun,costArray,iter] = update(obj,u,bc,costArray)
            i = 0; err = 1; costOld = costArray(end);

            LHS = obj.functional.computeHessian(u);
            RHS = obj.functional.computeGradient(u,bc);
            u0A = obj.computeDisplacement(LHS,RHS,u,bc);
            J0A = obj.functional.computeCost(u0A,bc)
            
            u0B = obj.updateDirichletValues(u,bc);
            J0B  = obj.functional.computeCost(u0B,bc)
            u = u0B;

            while (abs(err) > obj.tol) && (i < obj.maxIter)
                LHS = obj.functional.computeHessian(u);
                RHS = obj.functional.computeGradient(u,bc);
                u   = obj.computeDisplacement(LHS,RHS,u,bc);

                [err, cost] = obj.computeErrorCost(u,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;

                i = i+1;
                obj.monitor.printCost('iterU',i,cost,err);
                obj.monitor.update(length(costArray),{[],[cost],[],[]});
                obj.monitor.refresh(); 
            end
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
        end

        function u = computeDisplacement(obj,LHS,RHS,u,bc)
            [LHSr,RHSr] = obj.fullToReduced(LHS,RHS,bc);
            if ~isempty(LHSr)
                uV = reshape(u.fValues',[u.nDofs 1]);
                uNewV = uV;

                uInFree = uV(bc.free_dofs);
                uNewFree = obj.updateWithNewton(LHSr,RHSr,uInFree);
                uNewV(bc.free_dofs) = uNewFree;
                uNew = reshape(uNewV,[flip(size(u.fValues))])';
            else
                uNew = u.fValues;
            end
            u.setFValues(uNew);
        end

        function u = updateDirichletValues(~,u,bc)
            fixedDofs = bc.dirichlet_dofs;
            if ~isempty(fixedDofs)
                dirich = bc.dirichletFun;
                uV   = reshape(u.fValues',[u.nDofs 1]);
                dirichVec = reshape(dirich.fValues',[dirich.nDofs 1]);
                uV(fixedDofs) = dirichVec(fixedDofs);
                uV = reshape(uV,[flip(size(u.fValues))])';
            end
            u.setFValues(uV);
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