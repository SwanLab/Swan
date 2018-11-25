classdef LevelSetCircle < LevelSetSphereNdim
    
    methods (Access = public)
        
        function obj = LevelSetCircle(input)
            obj.fracRadius = 0.4;
            obj.compute(input);
        end
    end    
    
end


