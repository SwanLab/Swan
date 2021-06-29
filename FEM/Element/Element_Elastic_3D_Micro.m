classdef Element_Elastic_3D_Micro < Elastic3D & Element_Elastic_Micro
    
    methods (Access = public)
        function obj = Element_Elastic_3D_Micro(mesh,geometry,material,dof,problemData,interpU)
            obj.compute(mesh,geometry,material,dof,problemData,interpU);
        end
        
    end
   
end