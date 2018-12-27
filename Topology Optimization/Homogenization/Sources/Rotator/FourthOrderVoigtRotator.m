classdef FourthOrderVoigtRotator < Rotator
    
    properties
    end
    
    methods (Access = protected)
        
        function generateRotator(obj)
            a = obj.angle;
            d = obj.dir;
            rotator = obj.createStressRotator(a,d);
            obj.rotationMatrix = rotator.getRotationMatrix();
        end
        
        function computeRotation(obj,tensor)
            obj.createRotatedTensor(tensor);
            R = obj.rotationMatrix;
            C = tensor.getValue();
            Cr = R*C*R';
            obj.rotatedTensor.setValue(Cr);
        end
    end
    
    methods (Access = protected, Abstract,Static)
        createStressRotator(obj)
    end
    
end

