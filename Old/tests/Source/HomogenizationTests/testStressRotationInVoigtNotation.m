classdef testStressRotationInVoigtNotation < handle
    
    properties (Access = public)
        tol =  1e-14;
    end
    
    properties (Access = protected)
        direction
        stress
        stressVoigt
        rotatedStressByVoigt
        rotatedStress
    end
    
    properties (Access = private)
        angle
    end
    
    methods (Access = public)

        function error = computeError(obj)
            rotStre        = obj.rotatedStress.getValue();
            rotStreByVoigt = obj.rotatedStressByVoigt.getValue();
            error = norm(double(rotStre) - double(rotStreByVoigt));
        end

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
            obj.stress = Stress3DTensor;
            obj.stress.createRandomTensor();
            obj.stressVoigt = Tensor2VoigtConverter.convert(obj.stress);
        end
        
        function createRotatedStress(obj)
            rotS = obj.rotateStress();
            obj.rotatedStress = Tensor2VoigtConverter.convert(rotS);            
        end
        
        function rotS = rotateStress(obj)
            a  = obj.angle;
            d  = obj.direction;
            rotS = Rotator.rotate(obj.stress,a,d);
        end
        
    end
    
    methods (Access = private)
        
        function createAngle(obj)
            obj.angle = 2*pi*rand(1);
        end
        
        function createRotatedStressWithVoigtNotation(obj)
            theta = obj.angle;
            dir   = obj.direction;
            stre = obj.stressVoigt;
            rotStre = Rotator.rotate(stre,theta,dir);
            obj.rotatedStressByVoigt = rotStre;
        end
        
    end
    
    methods (Abstract,Access = protected)
        createDirection(obj)
    end
end

