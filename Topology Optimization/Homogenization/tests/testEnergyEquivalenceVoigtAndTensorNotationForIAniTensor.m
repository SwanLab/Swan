classdef testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor ...
        < testEnergyEquivalenceVoigtAndTensorNotation
    
    properties
    end
    
    methods (Access = public)
        
        function obj = testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor()
            obj@testEnergyEquivalenceVoigtAndTensorNotation();
        end 
        
    end
    
    methods        
        function generateFourthOrderTensor(obj)
            obj.Ch = fourthOrderTensor();
            obj.Ch.createRandomTensor();
            obj.Ch.computeTensorVoigt();
        end        
    end
    
end

