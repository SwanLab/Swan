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

        
        function [uOutNew,uOutVec] = computeDisplacement (obj,Ktan,costOld,u,bc)
            % res = costOld ???
            % gradient =costOld
            % falla pq crida a gradient() per algun motiu 
            

            u = obj.update(costOld,u); %    <------------- ERROR



            [err,~] = computeErrorCost(obj,u,u,bc,costOld);
            while(err>0 && ~obj.isTooSmall())
                u.recoverOld();
                obj.decreaseStepLength();
                u = obj.update(g,u);
                [err,~] = computeErrorCost(obj,u,u,bc,costOld);
            end
            obj.increaseStepLength(10);
            uOutNew = u;
            uOutVec = reshape(uIn.fValues', [uIn.nDofs 1]);
        end
        
      
        function u = update(obj,g,u)  
            y  = u.fValues;
                y = reshape(y.', [], 1);

            t  = obj.tau;
            
            y  = y - t*g;
                y = reshape(y, 2, []).';
            u.update(y);
        end

        function computeFirstStepLength(obj,g,x,f)
            xVal    = x.fun.fValues;
            obj.tau = min(f*sqrt(norm(g)/norm(xVal)),obj.tauMax);
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

        % function updateBounds(obj,ub,lb)
        %     if ~isnumeric(ub)
        %         obj.upperBound = ub.fValues;
        %     else
        %         obj.upperBound = ub;
        %     end
        %     if ~isnumeric(lb)
        %         obj.lowerBound = lb.fValues;
        %     else
        %         obj.lowerBound = lb;
        %     end
        % end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.tauMax     = cParams.tauMax;
            obj.tau        = cParams.tau;
            obj.functional = cParams.functional;
        end

        % function updateBoundsMultipliers(obj,x,y)
        %     t          = obj.tau;
        %     dyx        = y-x;
        %     dxy        = x-y;
        %     lUB        = zeros(size(x));
        %     lLB        = zeros(size(x));
        %     lUB(dyx>0) = dyx(dyx>0);
        %     lLB(dxy>0) = dxy(dxy>0);
        %     obj.boxConstraints.lUB    = lUB;
        %     obj.boxConstraints.lLB    = lLB;
        %     obj.boxConstraints.refTau = t;
        % end
        
        
        
        function [e, cost] = computeErrorCost(obj,u,phi,bc,costOld)
            cost = obj.functional.computeCost(u,phi,bc);
            e = cost - costOld;
        end
    end


    
end