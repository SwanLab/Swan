classdef testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor ...
        < testEnergyEquivalenceVoigtAndTensorNotation

    methods (Access = public)

        function obj = testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor()
            obj@testEnergyEquivalenceVoigtAndTensorNotation();
        end

    end

    methods

        function generateFourthOrderTensor(obj)
            obj.Ch = Stiffness3DTensor();
            obj.Ch.createRandomTensor();
        end

    end

end