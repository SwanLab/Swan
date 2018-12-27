classdef testInverseSymmetricFourthOrderTensor < testInverseFourthOrderTensor
    
    methods (Access = protected)
        function createRandomFourthOrderTensor(obj)
            obj.tensor = Stiffness3DTensor;
            obj.tensor.createRandomTensor();
        end
    end
    
    methods (Static, Access = protected)
        
        function Id = computeIdentityTensor(I,i,j,k,l)
            Id = 0.5*(I(i,k)*I(j,l) + I(i,l)*I(j,k));
        end
    end
    
end