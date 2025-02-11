classdef StressRotator < Rotator
    
    
    methods (Access = public)
        
        function obj = StressRotator(angle,dir)
            obj.compute(angle,dir)
        end
        
    end
    
    methods (Access = protected)
        function generateRotator(obj)
            a = obj.angle;
            d = obj.dir;
            rotator = VectorRotator(a,d);
            obj.rotationMatrix = rotator.getRotationMatrix();
        end

        function computeRotation(obj,tensor)
            obj.createRotatedTensor(tensor);
            R = obj.rotationMatrix;
            s = tensor.getValue();
            sR = R'*s*R;
            obj.rotatedTensor.setValue(sR);
        end
        
    end
    
end

