classdef LevelSetCircle < LevelSetSphereNdim
    

    
    methods (Access = public)
        
        function obj = LevelSetCircle(input)
            obj.fracRadius = 0.4;
            obj.compute(input);
        end
    end    
    
    methods (Access = protected)

        
        function computeDesignVariable(obj)
            phi = obj.levelSet;
             switch obj.optimizerName
                case {'SLERP','HAMILTON-JACOBI'}
                    obj.x = phi;
                otherwise
                    initial_holes = ceil(max(phi,0))>0;
                    obj.x = obj.ini_design_value*ones(obj.lsSize);
                    obj.x(initial_holes) = obj.hole_value;
             end                        
        end
        
    end
    
end


