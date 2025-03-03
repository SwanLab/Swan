classdef SamplePointsCreatorForOptimalExponentComputer < handle
    
    properties (Access = public)
       rhoV
       txiV
       phiV
    end
    
    methods (Access = public, Static)
       
        function obj = create(cParams)
            f = SamplePointsCreatorForOptimalExponentComputerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)       
        
        function compute(obj)
            obj.computeRhoTxiValues()
        end        
        
    end
    
    methods (Access = protected, Abstract)
        computeRhoTxiValues(obj)
    end
    
    
end