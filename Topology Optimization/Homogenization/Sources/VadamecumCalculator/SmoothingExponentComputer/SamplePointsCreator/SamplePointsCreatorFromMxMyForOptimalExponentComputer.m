classdef SamplePointsCreatorFromMxMyForOptimalExponentComputer ...
        < SamplePointsCreatorForOptimalExponentComputer
    
    properties (Access = private)
        phiV
        mxMax
        myMax
        mxV
        myV
        superellipse
    end
    
    methods (Access = public)
        
        function obj = SamplePointsCreatorFromMxMyForOptimalExponentComputer()
            obj.init()            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.mxMax = 0.99;
            obj.myMax = 0.99;
            nMx = 10;
            nMy = 10;
            nPhi = 10;
            obj.mxV = linspace(0.01,obj.mxMax,nMx);
            obj.myV = linspace(0.01,obj.myMax,nMy);            
            obj.phiV = linspace(pi/4,pi/4,nPhi);   
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