classdef FourthOrderVoigtPlaneStressRotator < Rotator
    
    properties (Access = protected)
        rotatedTensor
        rotationMatrix
    end
    
    methods (Access = public)
        
        function obj = FourthOrderVoigtPlaneStressRotator(angle,dir)
            obj.init(angle,dir)
            obj.generateRotator()
        end
        
    end
    
    methods (Access = private)
        
        function generateRotator(obj)
            a = obj.angle;
            d = obj.dir;
            rotator = StressVoigtPlaneStressRotator(a,d);
            obj.rotationMatrix = rotator.getRotationMatrix();
        end
        
    end
    
    methods (Access = protected)
        function computeRotation(obj,tensor)
            R = obj.rotationMatrix;
            C = tensor.getValue();
            obj.rotatedTensor = R*C*R';
        end
        
    end
    
end

