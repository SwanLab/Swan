classdef Element_Elastic_2D_Micro < Elastic2D & Element_Elastic_Micro
        
    methods (Access = public)
        
        function obj = Element_Elastic_2D_Micro(mesh,geometry,material,dof)
            obj.compute(mesh,geometry,material,dof);
        end        

    end
    
    
end
