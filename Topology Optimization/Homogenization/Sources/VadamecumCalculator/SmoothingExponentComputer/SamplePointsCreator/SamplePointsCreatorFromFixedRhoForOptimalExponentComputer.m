classdef SamplePointsCreatorFromFixedRhoForOptimalExponentComputer ...
        < SamplePointsCreatorForOptimalExponentComputer
    
    properties (Access = public)
        txiMax
        txiMin        
    end
    
    properties (Access = private)
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
            obj.psiV = cParams.psi;                        
            obj.mxMax = 0.99;
            obj.myMax = 0.99;       
            qMax = 32;
            sE = SuperEllipseParamsRelator;
            obj.txiMax = sE.txiFromMxRho(obj.mxMax,obj.rho,qMax);
            obj.txiMin = sE.txiFromMxRho(obj.myMax,obj.rho,qMax);
            obj.npoints = 20;
        end
        
    end
    
    methods (Access = protected)
        
        function computeRhoTxiValues(obj)
            obj.rhoV = obj.rho0*ones(obj.npoints,1);
            obj.txiV(:,1) = linspace(obj.txiMin,obj.txiMax,obj.npoints);    
        end
  
    end
    
end