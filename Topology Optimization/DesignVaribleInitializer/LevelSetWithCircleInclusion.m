classdef LevelSetWithCircleInclusion < LevelSetSphereNdim
   
    methods (Access = public)
        
        function obj = LevelSetWithCircleInclusion(input)
            obj.fracRadius = 0.4;
            obj.compute(input);
        end
    end
    
    
    methods (Access = protected)
        
        function computeDesignVariable(obj)            
            phi = obj.levelSet;
            phi = -phi;
            obj.x = obj.ini_design_value*ones(obj.lsSize);
            obj.x( phi > 0 ) = obj.hole_value;
        end
        
    end
    
end

