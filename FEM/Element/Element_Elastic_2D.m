classdef Element_Elastic_2D < Elastic2D & Element_Elastic

    
    methods (Access = public)
        function obj = Element_Elastic_2D(mesh,geometry,material,dof,problemData,interpU)
            obj.compute(mesh,geometry,material,dof,problemData,interpU);
            obj.dEps = obj.computedEps();            
        end
    end
    
end