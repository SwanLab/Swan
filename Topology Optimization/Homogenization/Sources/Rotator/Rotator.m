classdef Rotator < handle
    
    properties (Access = public)
        rotationMatrix
    end
    
    properties (Access = protected)
        rotatedTensor
    end
    
    properties (Access = protected)
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
        
        function compute(obj,angle,dir)
            obj.init(angle,dir)
            obj.generateRotator()
        end
        
        function init(obj,angle,direction)
             obj.angle = angle;
             obj.dir = direction;
        end
        
        function t = getRotatedTensor(obj)
            t = obj.rotatedTensor;
        end
        
        function createRotatedTensor(obj,tensor)
            obj.rotatedTensor = tensor.clone();
        end
        
    end
    
    methods (Abstract, Access = protected)
        computeRotation(obj,tensor)
        generateRotator(obj)
    end
    
end

