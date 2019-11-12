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
            qmax = 32;
            c = obj.cFunction(qmax);
            obj.txiMax = atan((obj.mxMax*c)/(1-obj.rho0));
            obj.txiMin = atan((1-obj.rho0)/(obj.myMax*c)); 
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