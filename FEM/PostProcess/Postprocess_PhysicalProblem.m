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
        

        
        
        function Print_results(obj,results)
            switch obj.ptype
                case 'ELASTIC'
                    obj.Print_results_mechanics(results);
                case 'Stokes'
                    obj.Print_results_fluids(results);
            end
            
  
            
            
        end
        
        function Print_results_mechanics(obj,results)
            obj.PrintVector(obj.displ_name,obj.displ_component,'Elastic Problem','Vector','OnNodes','',results.physicalVars.d_u);
            obj.PrintTensor(obj.stress_name,obj.stress_component,'Elastic Problem','Vector','OnGaussPoints',obj.gauss_points_name,results.physicalVars.stress);
            obj.PrintTensor(obj.strain_name,obj.strain_component,'Elastic Problem','Vector','OnGaussPoints',obj.gauss_points_name,results.physicalVars.strain);
        end
        
        function Print_results_fluids(obj,results)
            fprintf(obj.fid_res,'OnGroup "u" \n');
            obj.PrintVector(obj.velocity_name,obj.velocity_component,'Stokes problem','Vector','OnNodes','',results.physicalVars.u);
            fprintf(obj.fid_res,'end ongroup\n\n');
            fprintf(obj.fid_res,'OnGroup "p" \n');
            obj.PrintScalar(obj.pressure_name,obj.pressure_component,'Stokes problem','Scalar','OnNodes','',results.physicalVars.p);
            fprintf(obj.fid_res,'end ongroup\n');
        end

    end
    
end

