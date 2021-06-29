classdef PonderatedOptimalSuperEllipseComputer < handle
    
    properties (Access = public)
        qMean
        qDesv
    end
    
    properties (Access = private)
        gaussian
        phiV
        qV
        xiV
        rhoV
    end
    
    methods (Access = public)
        
        function obj = PonderatedOptimalSuperEllipseComputer(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            nPoints = size(obj.rhoV,1);
            obj.qMean  = zeros(nPoints,1);
            obj.qDesv  = zeros(nPoints,1);            
            for iPoint = 1:nPoints                
                q   = obj.qV(iPoint,:);
                phi = obj.phiV(iPoint,:);
                xi  = obj.xiV(iPoint,1);           
                [qM,qD] = obj.computeMeandAndDeviation(phi,xi,q);    
                obj.qMean(iPoint,1) = qM;
                obj.qDesv(iPoint,1) = qD;
            end
        end
        
    end    
    
    methods (Access = private)
        
        function init(obj,cParams)
            v = cParams.vademecum;
            obj.rhoV = v.rhoV;
            obj.xiV  = v.xiV;
            obj.phiV = v.phiV;
            obj.qV   = v.qV;
        end
        
    end
    
    methods (Access = private, Static)
        
        function [qM,qD] = computeMeandAndDeviation(phi,xi,q)
            s.phiV = phi;
            c = SuperEllipseMeanAndDesvExponentComputer(s);
            qM = c.computeMean(q,xi);
            qD = c.computeDesv(q,xi);
        end
       
    end

    
end