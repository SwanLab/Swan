classdef StressVoigtPlaneStressRotator < StressVoigtRotator
    
    
    methods (Access = public)
        
        function obj = StressVoigtPlaneStressRotator(angle,dir)
            obj.compute(angle,dir)
        end
        
    end
    
    methods (Access = protected)
        function generateRotator(obj)
            a = obj.angle;
            d = obj.dir;
            rotatorVoigt3D = StressVoigt3DRotator(a,d);
            rot = rotatorVoigt3D.getRotationMatrix();
            obj.rotationMatrix = obj.makeItPlaneStress(rot);
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

