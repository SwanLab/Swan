classdef ShFunc_Chomog_EnforceCh < ShFunc_Chomog
    
    properties (Access = public)
        ChTarget
    end
    
    properties (Access = private)
       C0
       weights
       difC
       difCp
    end
    
    properties (Access = protected)
        pNorm
    end
    
    methods (Access = public)
        
      function computeGradientValue(obj)
            obj.computeChDerivative();
            nElem  = obj.getnElem();
            nGaus  = obj.getnGaus();
            nStres = obj.getnStre();
            nComp  = obj.computeNcomp(nStres);
            grad = zeros(nElem,nGaus,nComp);
            dCh = permute(obj.dCh,[3,4,1,2]);
            for iStres = 1:nStres
                for jStres = 1:nStres
                    iv     = obj.vector2Voigt(iStres,jStres);
                    difCij = obj.difC(iStres,jStres);
                    C0ij   = obj.C0(iStres,jStres);
                    wij    = obj.weights(iv);
                    p      = obj.pNorm;
                    dChij   = squeeze(dCh(:,:,iStres,jStres));
                    dCostIv = (sum(obj.difCp))^(1/p-1)*wij*(difCij/C0ij)^(p-1)*(dChij/C0ij);
                    grad(:,:,iv) = dCostIv + grad(:,:,iv);
                end
            end
            obj.gradient = sum(grad,3);
        end
        
        function computeFunctionValue(obj)
            obj.difC  = obj.computeDifC();
            obj.difCp = obj.computeNormP(obj.difC);
            p = obj.pNorm;
            obj.value = sum(obj.difCp)^(1/p);
        end
        
        function v = getVariablesToPlot(obj)
            ChTargetP = obj.computeNormP(obj.ChTarget);
            p = obj.pNorm;
            normChTarget = sum(ChTargetP)^(1/p);
            v{1} = (obj.value*obj.value0)/normChTarget;
        end
        
    end
    
    methods (Access = protected)
        
        function difCp = computeNormP(obj,difC)
            nStres = obj.getnStre();
            nComp  = obj.computeNcomp(nStres);
            dCp = zeros(nComp,1);
            for iStres = 1:nStres
                for jStres = 1:nStres
                    iv = obj.vector2Voigt(iStres,jStres);
                    difCij = difC(iStres,jStres);
                    C0ij   = obj.C0(iStres,jStres);
                    wij    = obj.weights(iv);
                    p      = obj.pNorm;
                    costij = wij*(difCij/C0ij)^p;
                    dCp(iv) = costij;
                end
            end
            difCp = dCp;
        end
        
        function difC = computeDifC(obj)
            difC = (obj.Chomog - obj.ChTarget);
        end
        
        function computeChTarget(obj,cParams)
            obj.ChTarget = ChTargetFactory.create(cParams);
        end
        
        function computeC0(obj)
            C = ones(3,3);
            C(1,1) = obj.ChTarget(1,1);
%             C(2,2) = obj.ChTarget(2,2);
%             C(3,3) = obj.ChTarget(3,3);
             C(1,2) = 0.01;
             C(2,1) = 0.01;
            obj.C0 = C;
        end
        
        function computeWeights(obj)
            f = 1;%0000;
            obj.weights = [1,1,1,1,1,f]';
        end
        
    end
    
    methods (Access = private, Static)
        
        function n = computeNcomp(nStre)
           n = (nStre+1)*nStre/2;
        end
        
        function [iv] = vector2Voigt(iStre,jStre)
            T = [1 6 5; 6 2 4; 5 4 3];
            iv = T(iStre,jStre);
        end
        
    end
end