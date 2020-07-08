classdef SuperEllipseMeanAndDesvExponentComputer < handle
    
    properties (Access = private)
        phiV
    end
    
    methods (Access = public)
        
        function obj = SuperEllipseMeanAndDesvExponentComputer(cParams)
            obj.init(cParams)
        end
        
        function qmean = computeMean(obj,qOpt,xi)
            qL = obj.computeQmean(qOpt,xi);
            qR = obj.computeQmean(qOpt,pi - xi);
            qmean = 0.5*qL + 0.5*qR;
        end   
        
        function qDesv = computeDesv(obj,qOpt,xi)
            qL = obj.computeQdesv(qOpt,xi);
            qR = obj.computeQdesv(qOpt,pi - xi);
            qDesv = 0.5*qL + 0.5*qR;
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
        
        function qDesv = computeQdesv(obj,qOpt,xiV)
            qMean = obj.computeQmean(qOpt,xiV);
            P(:,1) = obj.obtainPvsPhi(obj.phiV,xiV);
            q(:,1) = qOpt;
            num = trapz(obj.phiV,P.*(q-qMean).^2);
            den = trapz(obj.phiV,P);
            qDesv(:,1) = num/den;            
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