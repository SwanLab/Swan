classdef PlottingVademecumSuperEllipseVsStress < handle
    
    properties (Access = private)
        phiV
        qOpt
        qMin
        qMax
        qMean
        xi
        rho
        compressedFileName
        vademecum
    end
    
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function obj = PlottingVademecumSuperEllipseVsStress()
            obj.init();
            obj.loadCompressedVademecum();
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.compressedFileName = 'OptimalSuperEllipses';
        end
        
        function loadCompressedVademecum(obj)
            d = load(obj.compressedFileName);
            obj.vademecum = d;
        end       
        
        function compute(obj)
            nMx  = size(obj.vademecum.q,1);
            nMy  = size(obj.vademecum.q,2);
            for iMx = 1:nMx
                for iMy = 1:nMy
                    obj.computeDataForPlotter(iMx,iMy)
                    obj.plot(iMx,iMy);
                end
            end
        end
        
        function computeDataForPlotter(obj,iMx,iMy)
            obj.computeXi(iMx,iMy);
            obj.computeRho(iMx,iMy);
            obj.computePhi(iMx,iMy);
            obj.computeQopt(iMx,iMy);
            obj.computeQminQmax();            
            obj.computeQmean();
        end
        
        function computeXi(obj,iMx,iMy)
            v = obj.vademecum;
            obj.xi = v.xi(iMx,iMy,1);
        end
                
        function computeRho(obj,iMx,iMy)
            v = obj.vademecum;
            obj.rho = v.rho(iMx,iMy,1);            
        end     
        
        function computePhi(obj,iMx,iMy)
            v = obj.vademecum;
            obj.phiV = squeeze(v.phi(iMx,iMy,:));
        end  
        
        function computeQopt(obj,iMx,iMy)
            v = obj.vademecum;
            obj.qOpt = squeeze(v.q(iMx,iMy,:));            
        end      
        
        function computeQminQmax(obj)
            s.mxMax = 0.99;
            s.myMax = 0.99;
            s.mxMin = 0.01;
            s.myMin = 0.01;
            s.xi =  obj.xi;
            s.rho = obj.rho;
            sC = SuperellipseExponentBoundsComputer(s);
            obj.qMin = sC.qMin;
            obj.qMax = sC.qMax;            
        end        
        
        function computeQmean(obj)
            s.phiV = obj.phiV;       
            qMeanC    = SuperEllipseMeanExponentComputer(s);
            obj.qMean = qMeanC.compute(obj.qOpt,obj.xi);            
        end
        
        function plot(obj,iMx,iMy)
            nMx  = size(obj.vademecum.q,1);            
            iter = (iMx-1)*nMx + iMy;
            s.phiV   = obj.phiV;
            s.qOpt   = obj.qOpt;
            s.qMin   = obj.qMin;
            s.qMax   = obj.qMax;
            s.qMean  = obj.qMean; 
            s.rho    = obj.rho;
            s.xi     = obj.xi;
            p = OptimalSuperEllipseExponentVsGaussianPlotter(s);
            p.plot(iter);
        end
        
%         
%         function isConsistent = checkConsistency(obj,qV,rhoV,xiV,mxV,myV,iMx,iMy,iphi)
%             xi  = xiV(iMx,iMy,iphi);
%             rho = rhoV(iMx,iMy,iphi);
%             q = qV(iMx,iMy,iphi);
%             mx = mxV(iMx);
%             my = myV(iMy);
%             errorMx = norm(mx - SuperEllipseParamsRelator.mx(xi,rho,q))/norm(mx)
%             errorMy = norm(my - SuperEllipseParamsRelator.my(xi,rho,q))/norm(my)
%             isMxEqual = errorMx < 1e-14;
%             isMyEqual = errorMy < 1e-14;
%             isConsistent = isMxEqual && isMyEqual;
%         end
        

        
    end
    
    
end