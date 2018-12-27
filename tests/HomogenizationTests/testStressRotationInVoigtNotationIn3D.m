classdef testStressRotationInVoigtNotationIn3D < testStressRotationInVoigtNotation
    
    properties (Access = private)

       
    end
    
    methods (Access = public)
        
        function obj = testStressRotationInVoigtNotationIn3D()
            obj.compute()
        end
        
    end

    methods (Access = protected)
        
        function createDirection(obj)
            obj.direction = Vector3D;
            dim = obj.direction.getTensorSize();
            obj.direction.setValue(rand(dim));
            obj.direction.normalize();
        end        

    end    
        
    
end