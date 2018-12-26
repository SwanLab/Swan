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
        
        function obj = ElasticityResultsPrinter(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results,iter)
            obj.init(fileID,fileName,nsteps,gaussDescriptor,etype,ptype,ngaus,ndim,posgp,results,iter)
            obj.print()
        end
    end
    
    methods (Access = protected)
        
        function printHeader(obj)
            obj.printGaussPointsHeader()
        end
        
        function printResults(obj)
            iS = obj.istep;
            gaussDescriptor = 'Guass up?';
            res = obj.results;
            VectorPrinter(obj.fileID,obj.displ_component,  res.physicalVars.d_u, obj.displ_name,iS,'OnNodes');
            TensorPrinter(obj.fileID,obj.stress_component, res.physicalVars.stress, obj.stress_name,iS,'OnGaussPoints',gaussDescriptor);
            TensorPrinter(obj.fileID,obj.strain_component, res.physicalVars.strain, obj.strain_name,iS,'OnGaussPoints',gaussDescriptor);
        end
        
    end
    
    methods (Access = private)
        
        function printGaussPointsHeader(obj)
            iD = obj.fileID;
            fprintf(iD,'GaussPoints "%s" Elemtype %s\n',obj.gaussDescriptor,obj.etype);
            fprintf(iD,'Number of Gauss Points: %.0f\n',obj.ngaus);
            fprintf(iD,'Nodes not included\n');
            fprintf(iD,'Natural Coordinates: given\n');
            for igaus = 1:obj.ngaus
                for idime = 1:obj.ndim
                    fprintf(iD,'%12.5d ',obj.posgp(igaus,idime));
                end
                fprintf(iD,'\n');
            end
            fprintf(iD,'End GaussPoints\n');
        end
        
    end
    
end

