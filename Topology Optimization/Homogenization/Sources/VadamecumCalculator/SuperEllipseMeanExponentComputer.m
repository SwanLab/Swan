classdef SuperEllipseMeanExponentComputer < handle
    
    properties (Access = private)
        phiV
    end
    
    methods (Access = public)
        
        function obj = SuperEllipseMeanExponentComputer(cParams)
            obj.init(cParams)
        end
        
        function qmean = compute(obj,qOpt,xi)
            qL = obj.computeQmean(qOpt,xi);
            qR = obj.computeQmean(qOpt,pi - xi);
            qmean = 0.5*qL + 0.5*qR;
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.phiV = cParams.phiV;
        end
        
        function qmean = computeQmean(obj,qOpt,xiV)
            P(:,1) = obj.obtainPvsPhi(obj.phiV,xiV);
            q(:,1) = qOpt;
            num = trapz(obj.phiV,P.*q);
            den = trapz(obj.phiV,P);
            qmean(:,1) = num/den;
        end        
   
    end
    
    methods (Access = public, Static)
        
        function P = obtainPvsPhi(phi,txi)            
            phiH = abs((phi) - (txi));
            phiS = abs((phi + pi) - (txi));
            phiT = [phiS,phiH];
            [~,ind] = min(phiT,[],2);            
            phiS = phi;
            phiS(ind == 1) = phi(ind == 1) + pi;
            phiS(ind == 2) = phi(ind == 2);
            gaussian = gaussianFunction();
            P = gaussian(txi,phiS);       
        end                
        
    end
    
end