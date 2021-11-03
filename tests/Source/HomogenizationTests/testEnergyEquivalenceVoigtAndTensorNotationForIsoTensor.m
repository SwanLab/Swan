classdef testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor  ...
         < testEnergyEquivalenceVoigtAndTensorNotation
    
    methods (Access = public)
        
        function obj = testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor()
            obj@testEnergyEquivalenceVoigtAndTensorNotation();
        end
        
    end

    methods 
        function generateFourthOrderTensor(obj)
             obj.Ch = IsotropicConstitutiveTensor(1,1/3);
        end
    end

end