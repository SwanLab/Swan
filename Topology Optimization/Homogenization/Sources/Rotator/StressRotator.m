classdef StressRotator < Rotator
    
    properties (Access = protected)
        rotatedTensor
        rotationMatrix
    end
    
    methods (Access = public)
        
        function obj = StressRotator(angle,dir)
            obj.init(angle,dir)
            obj.generateRotator()
        end
        
    end
    
    methods (Access = private)
        function generateRotator(obj)
            a = obj.angle;
            d = obj.dir;
            rotator = VectorRotator(a,d);
            obj.rotationMatrix = rotator.getRotationMatrix();
        end
        
    end
    
    methods (Access = protected)
        function computeRotation(obj,tensor)
            R = obj.rotationMatrix;
            s = tensor.getValue();
            obj.rotatedTensor = R'*s*R;
        end
        
    end
    
end

