classdef ChRotatorForFiberHomogenizer < handle
    
    methods (Access = public)
        
        function Ch = rotate(obj, dir,Ch)
            angle  = obj.computeRotationAngle(dir);
            dir    = obj.computeZdirection();
            C      = obj.makeChTensor(Ch);
            Ch = Rotator.rotate(C,angle,dir);
        end
        
    end
    
    methods (Access = private, Static)
        
         function angle = computeRotationAngle(dir)
            fiberDirection = dir.getValue;
            angle = -acos(dot(fiberDirection,[1 0 0]));
        end
        
        function C = makeChTensor(Ch)
            C = SymmetricFourthOrderPlaneStressVoigtTensor();
            C.setValue(Ch);
        end
        
        function dir = computeZdirection()
            d = [ 0 0 1];
            dir = Vector3D;
            dir.setValue(d);
        end
        
    end
end
