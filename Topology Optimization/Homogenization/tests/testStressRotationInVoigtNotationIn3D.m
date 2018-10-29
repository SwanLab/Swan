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
            obj.Direction = rand(3,1);
            obj.Direction = obj.Direction/norm(obj.Direction);
        end        

    end    
        
    
end