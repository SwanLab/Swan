classdef Postprocess_PhysicalProblem < Postprocess
    
    properties
        stress_name = 'Stress';
        stress_component = 'S';
        strain_name = 'Strain';
        strain_component = 'E';
        displ_name = 'Displacements';
        displ_component = 'U';
        velocity_name = 'Velocity';
        velocity_component = 'U'
        pressure_name = 'Pressure';
        pressure_component = 'p';
    end
    
    methods (Access = public)
        % Export Results
        

        
        
        function Print_results(obj,results,ifield,istep)
            switch obj.ptype
                case 'ELASTIC'
                    obj.Print_results_mechanics(results);
                case 'Stokes'
                    obj.Print_results_fluids(results,ifield,istep);
            end
            
  
            
            
        end
        
        function Print_results_mechanics(obj,results)
            obj.PrintVector(obj.displ_name,obj.displ_component,'Elastic Problem','Vector','OnNodes','',results.physicalVars.d_u);
            obj.PrintTensor(obj.stress_name,obj.stress_component,'Elastic Problem','Vector','OnGaussPoints',obj.gauss_points_name,results.physicalVars.stress);
            obj.PrintTensor(obj.strain_name,obj.strain_component,'Elastic Problem','Vector','OnGaussPoints',obj.gauss_points_name,results.physicalVars.strain);
        end
        
        function Print_results_fluids(obj,results,ifield,istep)
            if ifield == 1
                    obj.PrintVector(obj.velocity_name,obj.velocity_component,'Stokes problem','Vector','OnNodes','',results.physicalVars.u(:,istep),istep);
            else            
                    obj.PrintScalar(obj.pressure_name,obj.pressure_component,'Stokes problem','Scalar','OnNodes','',results.physicalVars.p(:,istep),istep);                
            end            
        end

    end
    
end

