classdef Postprocess_TopOpt_density < Postprocess_TopOpt
    properties
        density_name = 'density';
        density_name_component = 'Dens'
        
        density_name_reg = 'density_reg';
        density_name_component_reg = 'Dens_reg'
    end
    
    methods  (Access = public)        
        function Print_design_variable(obj,results)
            obj.Print_density(results);
        end
        
        function Print_density(obj,results)
            obj.PrintScalar(obj.density_name,obj.density_name_component,'Elastic Problem','Scalar','OnNodes','',results.design_variable,results.iter);
        end
        
        function Print_density_reg(obj,density)
            obj.PrintScalar(obj.density_name_reg,obj.density_name_component_reg,'Elastic Problem','Scalar','OnGaussPoints',obj.gauss_points_name,density);
        end
    end
end
