classdef SamplePointsCreatorFromMxMyForOptimalExponentComputer ...
        < SamplePointsCreatorForOptimalExponentComputer
    
    properties (Access = private)
        mxMax
        myMax
        mxV
        myV
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
            nPsi = 10;
            obj.mxV = linspace(0.01,obj.mxMax,nMx);
            obj.myV = linspace(0.01,obj.myMax,nMy);            
            obj.psiV = linspace(pi/4,pi/4,nPsi);            
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
                    c = obj.cFunction(q);
                    index = nmx*(imy - 1) + imx;
                    obj.rhoV(index) = 1 - c*mx*my;
                    obj.txiV(index) = atan(mx/my);
                end
            end            
        end
  
    end
    
    
end