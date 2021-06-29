classdef SamplePointsCreatorFromMxMyForOptimalExponentComputer ...
        < SamplePointsCreatorForOptimalExponentComputer
    
    properties (Access = private)        
        mxMax
        myMax
        mxMin
        myMin
        mxV
        myV
        superellipse
    end
    
    methods (Access = public)
        
        function obj = SamplePointsCreatorFromMxMyForOptimalExponentComputer(cParams)
            obj.init(cParams)            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            nMx  = cParams.nMx;
            nMy  = cParams.nMy;
            nPhi = cParams.nPhi;
            phiMin = cParams.phiMin;
            phiMax = cParams.phiMax;            
            obj.mxMax = 0.99;
            obj.myMax = 0.99;
            obj.mxMin = 0.01;
            obj.myMin = 0.01;
            obj.mxV = linspace(obj.mxMin,obj.mxMax,nMx);
            obj.myV = linspace(obj.myMin,obj.myMax,nMy);            
            obj.phiV = linspace(phiMin,phiMax,nPhi);   
            obj.superellipse = SuperEllipseParamsRelator();
        end
        
    end
    
    methods (Access = protected)
        
        function computeRhoTxiValues(obj)
            nmx = length(obj.mxV);
            nmy = length(obj.myV);
            q = 10^6;
            for imx = 1:nmx
                for imy = 1:nmy
                    mx = obj.mxV(imx);
                    my = obj.myV(imy);
                    index = nmx*(imy - 1) + imx;
                    obj.rhoV(index) = obj.superellipse.rho(mx,my,q);
                    obj.txiV(index) = obj.superellipse.txi(mx,my);
                end
            end            
        end
  
    end
       
    
end