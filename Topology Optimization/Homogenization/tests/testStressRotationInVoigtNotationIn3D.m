classdef testStressRotationInVoigtNotationIn3D < testStressRotationInVoigtNotation
    
    properties (Access = private)

       
    end
    
    methods (Access = public)
        
        function obj = testStressRotationInVoigtNotationIn3D()
            obj@testStressRotationInVoigtNotation()
        end
        
    end

    methods (Access = protected)
        
        function createDirection(obj)
            obj.Direction = rand(3,1);
            obj.Direction = obj.Direction/norm(obj.Direction);
        end        
        
        function setVoigtMatrix(obj,MatrixGenerator)
            obj.VoigtRotation = MatrixGenerator.VoigtMatrix;
        end
        
        function Tensor = obtainRotatedTensor(obj)
           Tensor = obj.RotatedStress.tensorVoigt; 
        end
        
        function Tensor = obtainRotatedTensorByVoigt(obj)
            Tensor = obj.RotatedStressByVoigt.tensorVoigt;            
        end
        
        function Stre = obtainStress(obj)
            Stre = obj.Stress.tensorVoigt;
        end
        
        function setRotatedStressByVoigt(obj,Tensor)
           obj.RotatedStressByVoigt.tensorVoigt = Tensor;
        end
        
    end    
        
    
end