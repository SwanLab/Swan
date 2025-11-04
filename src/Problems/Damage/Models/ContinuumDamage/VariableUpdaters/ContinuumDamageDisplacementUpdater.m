    classdef ContinuumDamageDisplacementUpdater < handle
        
        properties (Access = private)
            functional
            monitor
            mesh
    
            solverType
            solver
    
            tol
            maxIter
        end
    
        methods (Access = public)
    
            function obj = ContinuumDamageDisplacementUpdater(cParams)
                obj.init(cParams);
            end
    
            function [u,r,F,costArray,iter] = update(obj,u,r,bc,costArray)
                i = 0; err = 1; costOld = costArray(end);
                while (abs(err) > obj.tol) && (i < obj.maxIter)
                    tauEps = obj.functional.computeTauEpsilon(u);
                    r.update(tauEps);
    
                    [res]       = obj.functional.computeResidual(u,r,bc);
                    [Ktan,Ksec] = obj.functional.computeDerivativeResidual(u,r);
                    [uNew,uVec] = obj.solver.computeDisplacement(Ktan,res,u,bc);
                    u.setFValues(uNew);
    
                    [err, cost] = obj.computeErrorCost(u,r,bc,costOld);
                    costArray(end+1) = cost;
                    costOld = cost;
    
                    i = i+1;
                    obj.monitor.printCost('iter',i,cost,err);
                    obj.monitor.update(length(costArray),{[],[],[],[cost],[],[]});
                    obj.monitor.refresh();
                end
                r.updateRold();
                F = Ksec*uVec;
                iter = i;
            end
    
        end
    
        methods (Access = private)
    
            function init(obj,cParams)
                obj.functional = cParams.functional;
                obj.monitor    = cParams.monitor;
                obj.tol        = cParams.tolerance;
                obj.maxIter    = cParams.maxIter;
                obj.solverType = cParams.solverType;
                obj.setSolver();
            end
    
            % function [uOut,uOutVec] = computeDisplacement(obj,LHSfull, RHSfull,uIn,bc)
            %     [LHS,RHS] = fullToReduced(obj,LHSfull,RHSfull,bc);
            %     if ~isempty(LHS)
            %         uInVec = reshape(uIn.fValues',[uIn.nDofs 1]);
            %         uOutVec = uInVec;
            % 
            %         uInFree = uInVec(bc.free_dofs);
            %         uOutFree = obj.updateWithNewton(LHS,RHS,uInFree);
            %         uOutVec(bc.free_dofs) = uOutFree;
            %         uOut = reshape(uOutVec,[flip(size(uIn.fValues))])';
            %     else
            %         uOut = uIn.fValues;
            %         uOutVec = reshape(uIn.fValues',[uIn.nDofs 1]);
            %     end
            % end
    
            % function [LHS,RHS] = fullToReduced(~,LHS,RHS,bc)
            %     free_dofs = bc.free_dofs;
            %     LHS = LHS(free_dofs, free_dofs);
            %     RHS = RHS(free_dofs);
            % end
            % 
            % function xNew = updateWithNewton(~,LHS,RHS,x)
            %     deltaX = -LHS\RHS;
            %     xNew = x + deltaX;
            % end
    
            function [e, cost] = computeErrorCost(obj,u,r,bc,costOld)
                cost = obj.functional.computeEnergy(u,r,bc);
                e = cost - costOld;
            end
            
            function setSolver(obj)
                switch obj.solverType
                    case 'Newton'
                        obj.solver = Newton();
                    case 'AdaptiveGradient'
                        % s.ub            = 1;
                        % s.lb            = initialGuess.fun.fValues;%????????????
                        s.tauMax        = 1e10;
                        s.tau           = 100;%????????????
                        s.functional    =   obj.functional;
                        obj.solver = AdaptiveGradient(s);
                end
            end
    
        end
    
    end