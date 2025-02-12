classdef StressVoigtRotator < Rotator
    
    properties
    end
    
    methods (Access = protected)
        function computeRotation(obj,tensor)
            obj.createRotatedTensor(tensor);
            R = obj.rotationMatrix;
            s = tensor.getValue();
            sR = R*s;
            obj.rotatedTensor.setValue(sR);
        end
    end
    
end

