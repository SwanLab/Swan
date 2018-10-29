classdef testStressRotationInVoigtNotation < test
    
    properties (Access = protected)
        Direction
        Stress
        StressVoigt
        rotatedStressByVoigt
        rotatedStress
    end
    
    properties (Access = private)
        Angle
    end
    
    methods (Access = protected)
        
        function compute(obj)
            obj.createAngle()
            obj.createDirection()
            obj.createStress()
            obj.createRotatedStress()
            obj.createRotatedStressWithVoigtNotation()
        end
        
        function createStress(obj)
            obj.Stress = StressTensor();
            sv = Tensor2VoigtConverter.convert(obj.Stress);
            obj.StressVoigt = StressVoigtTensor;
            obj.StressVoigt.setValue(sv)
        end
        
        function createRotatedStress(obj)
            rotS = obj.rotateStress();
            obj.rotatedStress = obj.convertInVoigt(rotS);            
        end
        
        function rotS = rotateStress(obj)
            a  = obj.Angle;
            d  = obj.Direction;
            sR = Rotator.rotate(obj.Stress,a,d);
            rotS = StressTensor();
            rotS.setValue(sR);
        end
        
        function tensVoigt = convertInVoigt(obj,tens)
            tensVoigt = StressVoigtTensor();
            t = Tensor2VoigtConverter.convert(tens);
            tensVoigt.setValue(t)            
        end
        
        function hasPassed = hasPassed(obj)
            rotStre        = obj.rotatedStress.getValue();
            rotStreByVoigt = obj.rotatedStressByVoigt();
            hasPassed = norm(double(rotStre) - double(rotStreByVoigt)) < 1e-14;
        end
        
    end
    
    methods (Access = private)
        
        function createAngle(obj)
            obj.Angle = 2*pi*rand(1);
        end
        
        function createRotatedStressWithVoigtNotation(obj)
            theta = obj.Angle;
            dir   = obj.Direction;
            stre = obj.StressVoigt.getValue();
            
            stress = StressVoigtTensor();
            stress.setValue(stre);
            
            rotStre = Rotator.rotate(stress,theta,dir);
            obj.rotatedStressByVoigt = rotStre;            
        end
        
    end
    
    methods (Abstract,Access = protected)
        createDirection(obj)
    end
end

