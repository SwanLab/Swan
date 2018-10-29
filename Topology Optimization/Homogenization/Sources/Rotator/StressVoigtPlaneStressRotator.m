classdef StressVoigtPlaneStressRotator < Rotator
    
    properties (Access = protected)
        rotatedTensor
        rotationMatrix
    end
    
    methods (Access = public)
        
        function obj = StressVoigtPlaneStressRotator(angle,dir)
            obj.init(angle,dir)
            obj.generateRotator()
        end
        
    end
    
    methods (Access = private)
        function generateRotator(obj)
            a = obj.angle;
            d = obj.dir;
            rotatorVoigt3D = StressVoigtRotator(a,d);
            rot = rotatorVoigt3D.getRotationMatrix();
            obj.rotationMatrix = obj.makeItPlaneStress(rot);
        end
        
    end
    
    methods (Access = protected)        
        function computeRotation(obj,tensor)
            R = obj.rotationMatrix;
            s = tensor.getValue();
            obj.rotatedTensor = R*s;
        end        
    end
    
    methods (Access = private,Static)
        function r = makeItPlaneStress(rotation3D)
            psIndex = PlaneStressIndex();
            index(:,1) = psIndex.getInPlaneIndex();
            r = rotation3D(index,index);
        end
    end
    
    
end

