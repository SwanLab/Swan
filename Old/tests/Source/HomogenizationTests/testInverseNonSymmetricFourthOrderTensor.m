classdef testInverseNonSymmetricFourthOrderTensor < testInverseFourthOrderTensor

    methods (Access = protected)

        function createRandomFourthOrderTensor(obj)
            obj.tensor = FourthOrder3DTensor;
            obj.tensor.createRandomTensor();
        end

    end

    methods (Static, Access = protected)

        function Id = computeIdentityTensor(I,i,j,k,l)
           Id = I(i,k)*I(j,l);
        end
    end

end