classdef SamplePointsCreatorForOptimalExponentComputer < handle
    
    properties (Access = public)
       rhoV
       txiV
       psiV
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
    
    methods (Access = protected, Static)
        
        function c = cFunction(q)
            c = gamma(1 + 1/q)^2/gamma(1 + 2/q);
        end
        
    end        
    
end