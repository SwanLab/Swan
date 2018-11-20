classdef LevelSetWithSphereInclusion < LevelSetSphereNdim
    
    methods (Access = public)
        
        function obj = LevelSetSphere(input)
            obj.fracRadius = 1-1e-6;
            obj.compute(input);
        end
    end    
    
    methods (Access = protected)

        
        function computeDesignVariable(obj)
            phi = obj.levelSet;
            phi = -phi;
            obj.x = obj.ini_design_value*ones(obj.lsSize);
            obj.x(phi > 0) = obj.hole_value;
        end
        
    end
    
end


