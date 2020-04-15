classdef PonderatedOptimalSuperEllipseComputer < handle
    
    properties (Access = public)
        mxV
        myV
        txiV
        rhoV
        qMean
        qDesv
    end
    
    properties (Access = private)
        gaussian
        psiV
        qV
    end
    
    methods (Access = public)
        
        function obj = PonderatedOptimalSuperEllipseComputer()
            obj.readOptimalSuperEllipseData();
        end
        
        function compute(obj)
            obj.createGaussian();
            obj.computeQmean();
            obj.computeQdesv();
            obj.computeMxMy();            
        end
        
    end    
    
    methods (Access = private)
        
        function readOptimalSuperEllipseData(obj)
            qPath = 'Topology Optimization/Vademecums/';
            qFile = 'OptimalSuperEllipseExponentData';
            for i = 1:10
                fileName = [qPath,qFile,num2str(i),'.mat'];
                d = load(fileName);
                q(i,:) = d.x.q(i,:);
            end
            obj.rhoV = d.x.rho;
            obj.txiV = d.x.txi;
            obj.psiV = d.x.psi;
            obj.qV   = q;
        end
        
        function createGaussian(obj)
            obj.gaussian = gaussianFunction();
        end
        
        function computeQmean(obj)
            npoint = length(obj.rhoV);
            qmean = zeros(npoint,1);
            for ipoint = 1:npoint
                txi = obj.txiV(ipoint);
                psi = obj.psiV(:);
                q   = obj.qV(:,ipoint);
                P(:,1) = obj.gaussian(psi,txi);
                num = sum(P.*q);
                den = sum(P);
                qmean(ipoint,1) = num/den;
            end
            obj.qMean = qmean;
        end
        
        function computeQdesv(obj)
           npoint = length(obj.rhoV);
            qD = zeros(npoint,1);
            for ipoint = 1:npoint
                txi = obj.txiV(ipoint);
                psi = obj.psiV(:);
                q   = obj.qV(:,ipoint);
                qM  = obj.qMean(ipoint,1);
                int = (q - qM).^2;
                P(:,1) = obj.gaussian(psi,txi);
                num = sum(P.*int);
                den = sum(P);
                qD(ipoint,1) = num/den;
            end
            obj.qDesv = qD;            
        end
        
        function computeMxMy(obj)
            npoint = length(obj.rhoV);
            mx = zeros(npoint,1);
            my = zeros(npoint,1);
            for ipoint = 1:npoint
                rho = obj.rhoV(ipoint);
                txi = obj.txiV(ipoint);
                q = obj.qMean(ipoint,1);
                mx(ipoint) = obj.computeMx(txi,rho,q);
                my(ipoint) = obj.computeMy(txi,rho,q);
            end
            obj.mxV = mx;
            obj.myV = my;
        end
        
    end
    
    methods (Access = private, Static)
        
        function mx = computeMx(txi,rho,q)
            mx = SuperEllipseParamsRelator.mx(txi,rho,q);
        end
        
        function my = computeMy(txi,rho,q)
            my = SuperEllipseParamsRelator.my(txi,rho,q);            
        end
                
    end    
    
end