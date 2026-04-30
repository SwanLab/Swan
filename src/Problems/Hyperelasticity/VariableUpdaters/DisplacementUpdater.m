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

                fValuesArray(:,:,1) = u.fValues;%per comprobar

            while (abs(err) > obj.tol) && (i < obj.maxIter)
                LHS = obj.functional.computeHessian(u);
                RHS = obj.functional.computeGradient(u,bc);

                fValuesArray(:,:,i+2) = obj.computeDisplacement(LHS,RHS,u,bc); %per comprobar

                u.setFValues(obj.computeDisplacement(LHS,RHS,u,bc));

                [err, cost] = obj.computeErrorCost(u,bc,costOld);
                costArray(end+1) = cost;
                costOld = cost;

                i = i+1;
                % obj.monitor.printCost('iterU',i,cost,err);
                % obj.monitor.update(length(costArray),{[],[cost],[],[]});
                % obj.monitor.refresh(); 
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
        
            % % ---- Load / save LHS history ----
            % if isfile("LHSList.mat")
            %     S = load("LHSList.mat");
            %     LHSList = S.LHSList;
            % else
            %     LHSList = [];
            % end
            % i = size(LHSList,3) + 1;
            % LHSList(:,:,i) = LHS;
            % save("LHSList.mat", "LHSList");
            % % ---- Load / save RHS history ----
            % if isfile("RHSList.mat")
            %     S = load("RHSList.mat");
            %     RHSList = S.RHSList;
            % else
            %     RHSList = [];
            % end
            % i = size(RHSList,2) + 1;
            % RHSList(:,i) = RHS;
            % save("RHSList.mat", "RHSList");
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