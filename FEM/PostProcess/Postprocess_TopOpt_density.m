classdef Postprocess_TopOpt_density < Postprocess_TopOpt
    
    properties
        density_name = 'density';
        density_name_component = 'Dens'
        
        density_name_reg = 'density_reg';
        density_name_component_reg = 'Dens_reg'
        
    end
    
    methods  (Access = public)
        
         
         function obj = Postprocess_TopOpt_density()
             
             
         end
        
        function Print_results(obj,results)
            %Print Results
            obj.Print_results@Postprocess_PhysicalProblem(results)
            obj.Print_design_variable(results.design_variable)
            obj.Print_density_reg(results.design_variable_reg)
            
        end
        
    end
    
    methods (Access = protected)
        
        function Print_design_variable(obj,design_variable)
            obj.Print_density(design_variable);
        end
        
        function Print_density(obj,density)
            obj.PrintScalar(obj.density_name,obj.density_name_component,'Elastic Problem','Scalar','OnNodes','',density);
        end
        
        
         function Print_density_reg(obj,density)
            obj.PrintScalar(obj.density_name_reg,obj.density_name_component_reg,'Elastic Problem','Scalar','OnGaussPoints',obj.gauss_points_name,density);
        end
    end
    
    
    
end
