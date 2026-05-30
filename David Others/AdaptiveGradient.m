classdef AdaptiveGradient < handle
    
    properties (Access = public)
        tau
        boxConstraints
    end

    properties (Access = private)
        tauMax
        functional
    end
    
    methods (Access = public)
        function obj = AdaptiveGradient(cParams)
            obj.init(cParams);
        end

        
        function [uOutNew,uOutVec] = computeDisplacement (obj,LHSfull,RHSfull,u,bc,costOld,r)
            [~,RHS] = obj.fullToReduced(LHSfull, RHSfull, bc);
            uOld = u.fValues;
            u = obj.update(RHS,u,bc);

            [err,~] = computeErrorCost(obj,u,r,bc,costOld);
            while(err>0 && ~obj.isTooSmall())
                u.setFValues(uOld);
                obj.decreaseStepLength();
                u = obj.update(RHS,u,bc);
                [err,~] = computeErrorCost(obj,u,r,bc,costOld);
            end
            obj.increaseStepLength(10);
            uOutNew = u.fValues;
            uOutVec = reshape(u.fValues', [u.nDofs 1]);
        end
        
      
        function u = update(obj,RHS,u,bc)  
            uInVec = reshape(u.fValues', [u.nDofs 1]);
            uOutVec = uInVec;
            
            uInFree = uInVec(bc.free_dofs);
            uOutFree = obj.updateWithGradient(RHS,obj.tau,uInFree);

            uOutVec(bc.free_dofs) = uOutFree;
            uOut = reshape(uOutVec, [flip(size(u.fValues))])';

            u.setFValues(uOut);
        end

        
        function increaseStepLength(obj,f)
            obj.tau = min(f*obj.tau,obj.tauMax);
        end

        function decreaseStepLength(obj)
            obj.tau = obj.tau/2;
        end

        function is = isTooSmall(obj)
            is = obj.tau < 1e-10;
        end


    end

    methods (Access = private)
        function init(obj,cParams)
            obj.tauMax     = cParams.tauMax;
            obj.tau        = cParams.tau;
            obj.functional = cParams.functional;
        end

        function [e, cost] = computeErrorCost(obj,u,r,bc,costOld)
            cost = obj.functional.computeEnergy(u,r,bc);
            e = cost - costOld;
        end

        function [LHS, RHS] = fullToReduced(~,LHS, RHS, bc)
            free_dofs = bc.free_dofs;
            LHS = LHS(free_dofs, free_dofs);
            RHS = RHS(free_dofs);
        end
        
        function u = updateWithGradient(~,RHS,t,u) 
            u  = u - t*RHS;
        end


    end


    
end