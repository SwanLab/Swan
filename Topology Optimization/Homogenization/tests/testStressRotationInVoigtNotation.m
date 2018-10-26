classdef testStressRotationInVoigtNotation < test
    
    properties (Access = protected)
       Direction
       Rotation
       VoigtRotation
       Stress
       RotatedStressByVoigt
       RotatedStress
    end
    
    properties (Access = private)
       Angle
    end
    
    
    methods (Access = public)
        
        function obj = testStressRotationInVoigtNotation()
            obj.createAngle()
            obj.createDirection()
            obj.createRotationMatrix()
            obj.createVoigtRotationMatrix() 
            obj.createStress()
            obj.createRotatedStress()  
            obj.createRotatedStressWithVoigtNotation()
        end
        
    end
    
    methods (Access = private)
        
        function createAngle(obj)
             obj.Angle = 2*pi*rand(1);
        end
        
        function createRotationMatrix(obj)
           theta = obj.Angle;
           Dir   = obj.Direction;
           MatrixGenerator = RotationMatrixGenerator(theta,Dir);
           obj.Rotation = MatrixGenerator.Matrix();
        end
        
        function createVoigtRotationMatrix(obj)
           theta = obj.Angle;
           Dir   = obj.Direction;
           MatrixGenerator = VoigtRotationMatrixGenerator(theta,Dir);
           obj.setVoigtMatrix(MatrixGenerator)           
        end
        
        function createRotatedStressWithVoigtNotation(obj)
            Rot  = obj.VoigtRotation;
            Stre = obj.obtainStress();
            RotTensor = Rot*Stre;
            obj.RotatedStressByVoigt = StressTensor();            
            obj.setRotatedStressByVoigt(RotTensor)           
        end
        
    end
    
    methods (Abstract,Access = protected)
        createDirection(obj) 
        setVoigtMatrix(obj)
        obtainRotatedTensor(obj)
        obtainRotatedTensorByVoigt(obj)
        obtainStress(obj)
        setRotatedStressByVoigt(obj)
    end
     
    methods (Access = protected)
        
        function createStress(obj)
            obj.Stress = StressTensor();
            obj.Stress.transformTensor2Voigt();
        end
        
        function createRotatedStress(obj)
            Rot  = obj.Rotation;
            Stre = obj.Stress.tensor();
            obj.RotatedStress = StressTensor();
            obj.RotatedStress.tensor = Rot'*Stre*Rot;
            obj.RotatedStress.transformTensor2Voigt()
        end

        
        function hasPassed = hasPassed(obj)   
            RotStre        = obj.obtainRotatedTensor();
            RotStreByVoigt = obj.obtainRotatedTensorByVoigt();
            hasPassed = norm(double(RotStre) - double(RotStreByVoigt)) < 1e-14;
        end
        
    end
end

