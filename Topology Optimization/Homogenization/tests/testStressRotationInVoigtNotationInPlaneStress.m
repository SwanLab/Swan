classdef testStressRotationInVoigtNotationInPlaneStress < testStressRotationInVoigtNotation
    
    properties (Access = private)

       
    end
    
    methods (Access = public)
        
        function obj = testStressRotationInVoigtNotationInPlaneStress()
            obj@testStressRotationInVoigtNotation()
        end
        
    end
    
    methods (Access = protected)
        
        function createDirection(obj)
            obj.Direction = [0 0 1];
        end    
        
        function setVoigtMatrix(obj,MatrixGenerator)
            obj.VoigtRotation = MatrixGenerator.VoigtMatrixPlaneStress;
        end
        
        function createStress(obj)
            obj.createStress@testStressRotationInVoigtNotation()
            TensorPS = PlaneStressTransformer.transform(obj.Stress.tensorVoigt);
            obj.Stress.tensorVoigtInPlaneStress = TensorPS;
        end
        
        function createRotatedStress(obj)
            obj.createRotatedStress@testStressRotationInVoigtNotation()
            TensorPS = PlaneStressTransformer.transform(obj.RotatedStress.tensorVoigt);
            obj.RotatedStress.tensorVoigtInPlaneStress = TensorPS;
        end
            
        function Tensor = obtainRotatedTensor(obj)
           Tensor = obj.RotatedStress.tensorVoigtInPlaneStress; 
        end
        
        function Tensor = obtainRotatedTensorByVoigt(obj)
            Tensor = obj.RotatedStressByVoigt.tensorVoigtInPlaneStress;            
        end
        
        function Stre = obtainStress(obj)
            Stre = obj.Stress.tensorVoigtInPlaneStress;
        end            
        
        function setRotatedStressByVoigt(obj,Tensor)
           obj.RotatedStressByVoigt.tensorVoigtInPlaneStress = Tensor;
        end
            
    end
    
        
    
end