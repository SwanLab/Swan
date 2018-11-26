classdef Element_Elastic_2D < Elastic2D & Element_Elastic

    
    methods (Access = public)
        function obj = Element_Elastic_2D(mesh,geometry,material,dof)
            obj.compute(mesh,geometry,material,dof);
        end
    end
    
end