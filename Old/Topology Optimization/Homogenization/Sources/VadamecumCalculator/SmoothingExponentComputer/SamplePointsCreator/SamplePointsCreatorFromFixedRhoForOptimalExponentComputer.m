classdef SamplePointsCreatorFromFixedRhoForOptimalExponentComputer ...
        < SamplePointsCreatorForOptimalExponentComputer

    properties (Access = private)
        txiMax
        txiMin           
        rho0
        mxV
        myV
        mxMax
        myMax
        npoints
    end
    
    methods (Access = public)
        
        function obj = SamplePointsCreatorFromFixedRhoForOptimalExponentComputer(cParams)
            obj.init(cParams)            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.rho0 = cParams.rho0;
            obj.phiV = cParams.phi;                        
            obj.mxMax = 0.99;
            obj.myMax = 0.99;       
            qMax = 32;
            sE = SuperEllipseParamsRelator;
            obj.txiMax = sE.txiFromMxAndRho(obj.mxMax,obj.rho0,qMax);
            obj.txiMin = sE.txiFromMyAndRho(obj.myMax,obj.rho0,qMax);
            obj.npoints = 4;
        end
        
    end
    
    methods (Access = protected)
        
        function computeRhoTxiValues(obj)
            obj.rhoV = obj.rho0*ones(obj.npoints,1);
            obj.txiV(:,1) = linspace(obj.txiMin,obj.txiMax,obj.npoints);    
        end
  
    end
    
end