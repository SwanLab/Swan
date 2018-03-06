classdef Element_Hyperelastic_3D < Element_Hyperelastic
   
    properties
    end
    
    methods
        function obj = Element_Hyperelastic_3D(geometry)
            obj@Element_Hyperelastic(geometry); % when the child is created before the father
        end
        function obj = updateCoord(obj,u)
            % Update coordinates
            coord0 = obj.coord;
            coord  = reshape(coord0',[],1);
            coord(obj.dof.free) = coord(obj.dof.free) + u;
            coord0 = reshape(coord,3,[])';
            obj.coord = coord0;            
        end
    end
    
end