classdef Rotator < handle
    
    properties (Abstract, Access = protected)
        rotatedTensor
        rotationMatrix
    end
    
    properties (Access = protected)
       tensor
       angle
       dir
    end
    
    methods (Access = public,Static)
        
        function rotatedTensor = rotate(tensor,angle,direction)
            factory = RotatorFactory();
            rotator = factory.create(tensor,angle,direction);
            rotator.computeRotation(tensor)
            rotatedTensor = rotator.getRotatedTensor();
        end  
        
    end
    
    methods (Access = public)
        function r = getRotationMatrix(obj)
            r = obj.rotationMatrix;
        end
    end
    
    methods (Access = protected)
        
        function init(obj,angle,direction)
             obj.angle = angle;
             obj.dir = direction;
        end        
        
        function t = getRotatedTensor(obj)
            t = obj.rotatedTensor;
        end                
        
    end
    
    methods (Abstract, Access = protected)
        computeRotation(obj,tensor)
    end
    
end

