classdef Element_Elastic_2D_Micro < Element_Elastic_2D & Element_Elastic_Micro
    
    
    methods (Access = public)
        function obj = Element_Elastic_2D_Micro(mesh,geometry,material,dof)
            obj = obj@Element_Elastic_2D(mesh,geometry,material,dof);
        end
        
        function variables = computeVars(obj,uL)
            variables = computeVars@Element_Elastic_2D(obj,uL);
            variables = obj.computeStressStrainAndCh(variables);
        end
    end
    
    methods(Access = protected)
        
        function FextVolumetric = computeVolumetricFext(obj)
            FextVolumetric = computeVolumetricFext@Element_Elastic(obj);
            F_def = obj.computeStrainRHS(obj.vstrain);
            FextVolumetric = FextVolumetric + F_def;
        end
    end
    
    
end
