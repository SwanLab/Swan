classdef DensityGaussResultsPrinter < ResultsPrinter ...
    
    properties
        fieldName = 'RegularizedDensity';
        simulationCase = 'DensityGauss';
        headPrinter = GaussHeadPrinter;
    end
    
    methods (Access = public)
        
        function obj = DensityGaussResultsPrinter()
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
            dens = obj.fields.density; 
            iS = obj.istep;
            gaussDescriptor = 'Guass up?';
            dS = obj.createScalarDataBase(obj.fileID,dens, obj.fieldName,iS,'OnGaussPoints',gaussDescriptor);
            ScalarGaussPrinter(dS);            
        end
    
    end
    
    methods (Access = private)
        
        function d = createScalarDataBase(obj,fileID,fieldValues,fieldName,istep,fieldPosition,gaussDescriptor)
            d.fileID = fileID;
            d.fieldValues = fieldValues;
            d.fieldName = fieldName;
            d.istep = istep;
            d.fieldPosition = fieldPosition;
            d.gaussDescriptor = gaussDescriptor;
            d.simulationCase = obj.simulationCase;

        end        
        
    end
    
    
end