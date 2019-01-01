classdef ElasticityResultsPrinter < ResultsPrinter 
    properties (Access = private)
        stress_name = 'Stress';
        strain_name = 'Strain';
        displ_name  = 'Displacements';
        stress_component = 'S';
        strain_component = 'E';
        displ_component  = 'U';
        simulationCase = 'ElasticityResults';
        dSig
        dStr
        dV
        headPrinter = GaussHeadPrinter;
    end
    
    methods (Access = public)
        
        function obj = ElasticityResultsPrinter()
        end
    end
    
    methods (Access = protected)
        
        function printHeader(obj)
            d.fileID = obj.fileID;
            d.gaussDescriptor = obj.gaussDescriptor;
            d.etype = obj.etype;
            d.ngaus = obj.ngaus;
            d.ndim  = obj.ndim;
            d.posgp = obj.posgp;
            obj.headPrinter.print(d);
        end
        
        function printResults(obj)
            obj.createDataBases();
            VectorPrinter(obj.dV);
            TensorPrinter(obj.dSig);
            TensorPrinter(obj.dStr);
        end                
        
    end
    
    methods (Access = private)
        
        function createDataBases(obj)
            iS = obj.istep;
            gaussDescriptor = 'Guass up?';
            f = obj.fields;            
            obj.dV = obj.createVectorDataBase(obj.fileID,obj.displ_component,f.d_u,obj.displ_name,iS,'OnNodes');
            obj.dSig = obj.createTensorDataBase(obj.fileID,obj.stress_component, f.stress, obj.stress_name,iS,'OnGaussPoints',gaussDescriptor);
            obj.dStr = obj.createTensorDataBase(obj.fileID,obj.strain_component, f.strain, obj.strain_name,iS,'OnGaussPoints',gaussDescriptor);
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
            d.simulationCase = obj.simulationCase;
        end
        
        function d = createTensorDataBase(obj,fileID,fieldComponentName,fieldValues,fieldName,istep,fieldPosition,gaussDescriptor)
           d = obj.createVectorDataBase(fileID,fieldComponentName,fieldValues,fieldName,istep,fieldPosition);
           d.gaussDescriptor = gaussDescriptor;
        end     

        
    end
    
end

