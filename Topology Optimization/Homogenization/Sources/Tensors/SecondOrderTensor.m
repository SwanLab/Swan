classdef SecondOrderTensor < Tensor
    
    properties (Access = public)
        tensorVoigt
        tensorVoigtInPlaneStress
    end
    
    properties (Abstract)
        VoigtFactor
    end
    
    methods
        
        function obj = SecondOrderTensor()
            obj.generateRandomTensor();
            obj.transformTensor2Voigt();
        end
        
        function transformTensor2Voigt(obj)
            obj.tensorVoigt = Tensor2VoigtConverter.convert(obj);
        end
        

        function transformVoigt2Tensor(obj)
            obj.tensor = zeros(3,3);
            transformer = TensorVoigtIndexTransformer();
            for iVoigt = 1:length(obj.tensorVoigt)
                    [iv,jv] = transformer.voigt2tensor(iVoigt);
                    obj.computeVoigtFactor(iVoigt);
                    obj.tensor(iv,jv) = 1/obj.VoigtFactor*obj.tensorVoigt(iVoigt);
                    obj.tensor(jv,iv) = 1/obj.VoigtFactor*obj.tensorVoigt(iVoigt);                    
            end
        end
    end
    

    methods (Access = private)
        
        function generateRandomTensor(obj)
            obj.tensor = rand(3,3);
            obj.tensor = obj.symmetrize(obj.tensor);
        end
        
        function createTensorVoigt(obj)
            isNotSymbolic = isUnit(obj.tensor);
            if isNotSymbolic
               obj.tensorVoigt = zeros(6,1);
            else
               obj.tensorVoigt = sym(zeros(6,1));                    
            end            
        end
        
    end
    
    
    methods (Abstract)
        computeVoigtFactor(obj,a)
    end
 
    methods (Static)
        
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

