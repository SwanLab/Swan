classdef OppositeOrientationInterpolator < handle
    
    
    properties (Access = private)
       mesh
       orientationVector
    end
    
    methods (Access = public)
        
        function obj = OppositeOrientationInterpolator(cParams)
            obj.init(cParams);
        end
        
        function c = compute(obj)
            c = obj.computeSymmetricCondition();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh              = cParams.mesh;
            obj.orientationVector = cParams.orientationVector;
        end
        

    end

    
end