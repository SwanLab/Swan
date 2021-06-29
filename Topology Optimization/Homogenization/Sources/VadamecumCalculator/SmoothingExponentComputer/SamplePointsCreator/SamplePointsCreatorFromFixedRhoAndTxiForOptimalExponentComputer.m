classdef SamplePointsCreatorFromFixedRhoAndTxiForOptimalExponentComputer ...
        < SamplePointsCreatorForOptimalExponentComputer
    
    methods (Access = public)
        
        function obj = SamplePointsCreatorFromFixedRhoAndTxiForOptimalExponentComputer(cParams)
            obj.init(cParams)            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.rhoV = cParams.rho0;
            obj.phiV = cParams.phi;                        
            obj.txiV = cParams.txi;
        end
        
    end
    
    methods (Access = protected)
        
        function computeRhoTxiValues(obj)
        end
  
    end
    
end