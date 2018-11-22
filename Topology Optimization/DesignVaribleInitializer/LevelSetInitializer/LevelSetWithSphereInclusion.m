classdef LevelSetWithSphereInclusion < LevelSetSphereNdim
    
    methods (Access = public)
        
        function obj = LevelSetSphere(input)
            obj.fracRadius = 1-1e-6;
            obj.compute(input);
        end
    end    
    
    
end


