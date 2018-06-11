classdef Postprocess_PhysicalProblem < Postprocess
    
    properties
        stress_name = 'Stress';
        stress_component = 'S';
        strain_name = 'Strain';
        strain_component = 'E';
        displ_name = 'Displacements';
        displ_component = 'U';
    end
    
    methods (Access = public)
        % Export Results
        
        function obj = Postprocess_PhysicalProblem()
            
        end
        
        
        function Print_results(obj,results,iter)
            
            
            %% Print Results
            obj.PrintVector(obj.displ_name,obj.displ_component,'Elastic Problem','Vector','OnNodes','',results.physicalVars.d_u(:,iter));
            obj.PrintTensor(obj.stress_name,obj.stress_component,'Elastic Problem','Vector','OnGaussPoints',obj.gauss_points_name,results.physicalVars.stress(iter,:,:));
            obj.PrintTensor(obj.strain_name,obj.strain_component,'Elastic Problem','Vector','OnGaussPoints',obj.gauss_points_name,results.physicalVars.strain(iter,:,:));
        end
        

    end
    
end

