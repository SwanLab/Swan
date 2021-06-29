classdef ShFunc_Chomog_EnforceCh < ShFunc_Chomog
    
    properties (Access = public)
        ChTarget
    end
    
    
    methods (Access = protected)
        
        function filterGradient(obj)
            g = obj.gradient;
            gf = zeros(size(obj.Msmooth,1),obj.nVariables);
            for ivar = 1:obj.nVariables
                gs = g(:,:,ivar);
                gf(:,ivar) = obj.filter.getP1fromP0(gs);
            end
            %gf = obj.Msmooth*gf;
            g = gf(:);
            obj.gradient = g;
        end        
        
        function obj = passFilter(obj)
            mass=obj.Msmooth;
            g = obj.gradient;
            obj.gradient = zeros(size(mass,1),size(g,2));
            for t=1:size(obj.gradient,2)
                gradient=obj.filter.getP1fromP0(g(:,:,t));
                obj.gradient(:,:,t) = mass*gradient;
            end
            
        end
        
        
        %
        %         function computeCost(obj)
        %
        %         end
        
        function computeCCstar(obj)
            obj.compute_Chomog_Derivatives();
            
            %nChTarget = 1;%obj.computeL2Norm(obj.ChTarget);
            obj.Chomog
            obj.ChTarget
            incC = (obj.Chomog - obj.ChTarget)
            
            f = 1;%0000;
            weights = [1,1,1,1,1,f]';
            
            %Gradient
            neq = 6;
            nelem = obj.physicalProblem.element.nelem;
            ngaus = size(obj.tstrain,2);
            nStres = obj.physicalProblem.element.getNstre;
            cost = zeros(neq,1);            
            grad = zeros(nelem,ngaus,neq);            
            Cdiv = ones(3,3);
            Cdiv(1,1) = obj.ChTarget(1,1);
            Cdiv(2,2) = obj.ChTarget(2,2);
            Cdiv(3,3) = obj.ChTarget(3,3);
            
            for iStres = 1:nStres
                for jStres = 1:nStres
                    iv = obj.vector2Voigt(iStres,jStres);
                    costIv = weights(iv)*(incC(iStres,jStres)/Cdiv(iStres,jStres))^2;
                    DCij = squeeze(obj.Chomog_Derivatives(iStres,jStres,:,:));
                    dCostIv = 2*weights(iv)*incC(iStres,jStres)*DCij/Cdiv(iStres,jStres)^2;                    
                    cost(iv)         = costIv;
                    grad(:,:,iv) = dCostIv + grad(:,:,iv);                    
                end
            end
            obj.value = sum(cost);
            obj.gradient = sum(grad,3);
        end
        
        function computeChTarget(obj,cParams)
            obj.ChTarget = ChTargetFactory.create(cParams);
        end
        
    end
    
    methods (Access = private, Static)
        
        function [iv] = vector2Voigt(iStre,jStre)
            T = [1 6 5; 6 2 4; 5 4 3];
            iv = T(iStre,jStre);
        end
        
    end
end