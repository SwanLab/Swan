classdef SecondOrder3DTensor < AbstractTensor ...
                               & SecondOrderDescriptor ...
                               & TensorRepresentation ...
                               & Elasticity3dDescriptor
    
    methods (Access = public)
        
        function obj = SecondOrder3DTensor()
        end
        
        function createRandomTensor(obj)
            obj.createRandomTensor@AbstractTensor()
            tSym = obj.symmetrize(obj.getValue());
            obj.setValue(tSym)
        end
        
    end
    
    methods (Access = protected)

        function loadTensorSize(obj)
            obj.tensorSize = [3,3];
        end
    end
    
    methods (Static, Access = public)
        
         function isSymmetric = isSymmetric(A)
            Asym = SecondOrderTensor.symmetrize(A);
            isSymmetric = norm(Asym(:)) - norm(A(:)) < 1e-14;
        end

        function Asym = symmetrize(A)
            Asym = 0.5*(A+A');
        end
        
        function Asym = symmetrizeWithUpperDiagonal(A)
            Asym = A;
            Asym(2,1) = A(1,2);
            Asym(3,2) = A(2,3);
            Asym(3,1) = A(1,3);
        end
        
    end
    
end

