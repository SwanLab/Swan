classdef ElasticityResultsPrinter < ResultsPrinter
    
    properties (Access = private)%(GetAccess = protected, SetAccess = private)
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
        
        function obj = ElasticityResultsPrinter(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results)
            obj.init(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results)
            obj.print()
        end
    end
    
    methods (Access = protected)
        
        function printResults(obj,ifield,istep)
            switch obj.ptype
                case 'ELASTIC'
                    obj.Print_results_mechanics(obj.results,1);
                case 'Stokes'
                    obj.Print_results_fluids(obj.results,ifield,istep);
            end
        end
        
        function Print_results_mechanics(obj,results,istep)
            gaussDescriptor = 'Guass up?';
            VectorPrinter(obj.fileID,obj.displ_component, results.physicalVars.d_u, obj.displ_name,istep,'OnNodes');
            TensorPrinter(obj.fileID,obj.stress_component, results.physicalVars.stress, obj.stress_name,istep,'OnGaussPoints',gaussDescriptor);
            TensorPrinter(obj.fileID,obj.strain_component, results.physicalVars.strain, obj.strain_name,istep,'OnGaussPoints',gaussDescriptor);
            %obj.PrintVector(obj.displ_name,obj.displ_component,'Elastic Problem','Vector','OnNodes','',results.physicalVars.d_u,istep);
           % obj.PrintTensor(obj.stress_name,obj.stress_component,'Elastic Problem','Vector','OnGaussPoints',obj.gauss_points_name,results.physicalVars.stress,istep);
            %obj.PrintTensor(obj.strain_name,obj.strain_component,'Elastic Problem','Vector','OnGaussPoints',obj.gauss_points_name,results.physicalVars.strain,istep);
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

