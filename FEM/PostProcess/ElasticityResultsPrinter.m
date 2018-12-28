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
        
        function obj = ElasticityResultsPrinter()
        end
    end
    
    methods (Access = protected)
        
        function printHeader(obj)
            obj.printGaussPointsHeader()
        end
        
        function printResults(obj)
            iS = obj.istep;
            gaussDescriptor = 'Guass up?';
            f = obj.fields;
            dV = obj.createVectorDataBase(obj.fileID,obj.displ_component,f.d_u,obj.displ_name,iS,'OnNodes');
            dSig = obj.createTensorDataBase(obj.fileID,obj.stress_component, f.stress, obj.stress_name,iS,'OnGaussPoints',gaussDescriptor);
            dStr = obj.createTensorDataBase(obj.fileID,obj.strain_component, f.strain, obj.strain_name,iS,'OnGaussPoints',gaussDescriptor);
            VectorPrinter(dV);
            TensorPrinter(dSig);
            TensorPrinter(dStr);
        end
        
        
        
    end
    
    methods (Access = private)
        
        function d = createVectorDataBase(obj,fileID,fieldComponentName,fieldValues,fieldName,istep,fieldPosition)
            d.fileID = fileID;
            d.fieldComponentName = fieldComponentName;
            d.fieldValues = fieldValues;
            d.fieldName = fieldName;
            d.istep = istep;
            d.fieldPosition = fieldPosition;
        end
        
        function d = createTensorDataBase(obj,fileID,fieldComponentName,fieldValues,fieldName,istep,fieldPosition,gaussDescriptor)
            d.fileID = fileID;
            d.fieldComponentName = fieldComponentName;
            d.fieldValues = fieldValues;
            d.fieldName = fieldName;
            d.istep = istep;
            d.fieldPosition = fieldPosition;
            d.gaussDescriptor = gaussDescriptor;
        end

        
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

