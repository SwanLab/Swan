classdef AmplificatorTensorFromVademecum < VariableFromVademecum
    
    
    properties (Access = public)
       P2t  
    end
    
    properties (Access = protected)
        fieldName = 'Ptensor';
    end
    
    properties (Access = private)
        monomials
    end
    
    methods (Access = public)
        
        function [P,dP,P2,dP2] = compute(obj,x,iP)
            obj.computeParamsInfo(x);
            obj.setValuesToInterpolator(x);
            [P,dP,P2,dP2] = obj.computeValues(iP);
        end
        
    end
    
    methods (Access = protected)
        
        function obtainValues(obj)
            iPnorm = 1;
            obj.monomials = obj.vadVariables.monomials{iPnorm};
            var = obj.vadVariables.variables;
            mxV = obj.vadVariables.domVariables.mxV;
            myV = obj.vadVariables.domVariables.mxV;
            for imx = 1:length(mxV)
                for imy = 1:length(myV)
                    P = var{imx,imy}.(obj.fieldName);
                    Ptensor = obj.transformVector2tensor(P{iPnorm});
                    v(:,:,imx,imy) = Ptensor;
                end
            end
            obj.values = v;
            
           for imx = 1:length(mxV)
                for imy = 1:length(myV)
                    P = var{imx,imy}.(obj.fieldName);
                    for iPnorm = 1:numel(P)
                    obj.P2t{iPnorm}(imx,imy,:) = P{iPnorm};
                    end
                end
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function [P2,dP2,Pp,dPp] = computeValues(obj,pNorm)
            iP = obj.computePNormIndex(pNorm);
            nstre = size(obj.values,1);
            P2  = zeros(nstre,nstre,obj.nPoints);
            dP2 = zeros(nstre,nstre,obj.nPoints,obj.nParams);
            for i = 1:nstre
                for j = 1:nstre
                    pv = squeeze(obj.values(i,j,:,:));
                    [p,dp] = obj.interpolator.interpolate(pv);
                    P2(i,j,:) = p;
                    dP2(i,j,:,1) = dp(:,1);
                    dP2(i,j,:,2) = dp(:,2);
                end
            end
            ptt = obj.makeSmall(obj.P2t{iP},pNorm);
            [Pp0,dPp] = obj.interpolator.interpolate(ptt);
            Pp  = obj.makeLarge(Pp0,pNorm);
            dPp = obj.makeLarge2(dPp,Pp0,pNorm);
        end
        
        function y = symlog(obj,x)
            y = sign(x).*log(abs(x)+1);
        end
        
        function y = invsymlog(obj,x)
            y = sign(x).*exp(abs(x)-1);
        end
        
        function pt = makeSmall(obj,pt,pNorm)
          %  pt = pt;
          %  pt = sign(pt).*(abs(pt).^(1/(1*pNorm)));
            pt = obj.symlog(pt);
        end
        
        function Pp = makeLarge(obj,Pp,pNorm)
            Pp = obj.invsymlog(Pp);
           % Pp = sign(Pp).*(exp(1).^abs(Pp));
          %  Pp = sign(Pp).*(abs(Pp).^(1*pNorm));
          %  Pp = Pp;
        end
        
        function dPp = makeLarge2(obj,dPp,Pp,pNorm)
           % dPp = obj.invsymlog(dPp);
            dPp(:,:,1) = dPp(:,:,1).*exp(abs(Pp));
            dPp(:,:,2) = dPp(:,:,2).*exp(abs(Pp));
          %  dPp = dPp;
        end
        
        
        function iP = computePNormIndex(obj,p)
           pC = 2:2:32;
           isIndex = p == pC;
           iP = find(isIndex);
        end

        function P = transformVector2tensor(obj,Pv)
            nstre = 3;
            P = zeros(nstre,nstre);
            for i = 1:nstre
                for j = 1:nstre
                    k = obj.voigt2TensorWithMonomials(i,j);
                    if i~= j
                        f = 0.5;
                    else
                        f = 1;
                    end
                    P(i,j) = f*Pv(k);
                end
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function tij = voigt2TensorWithMonomials(i,j)
            T = [6 1 2;
                 1 5 3
                 2 3 4];
            tij = T(i,j);
        end
        
    end
    
end
