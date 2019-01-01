classdef DensityResultsPrinter < ResultsPrinter
    
    properties (Access = private)
        fieldName      = 'Density';
        simulationCase = 'NodalDensity';
        headPrinter = NoGaussHeadPrinter;
    end
    
    methods (Access = public)
        
        function obj = DensityResultsPrinter()
        end
    end
    
    
    methods (Access = protected)
        
        function printHeader(obj)
            obj.headPrinter.print(obj.fileID);
        end
        
        function printResults(obj)
            dens = obj.fields; 
            iS = obj.istep;
            dS = obj.createScalarDataBase(obj.fileID,dens,obj.fieldName,iS,'OnNodes');   
            ScalarPrinter(dS);
        end
    
    end
    
    methods (Access = private)
        
        function d = createScalarDataBase(obj,fileID,fieldValues,fieldName, istep, fieldPosition)
            d.fileID = fileID;
            d.fieldValues = fieldValues;
            d.fieldName = fieldName;
            d.istep = istep;
            d.fieldPosition = fieldPosition;
            d.simulationCase = obj.simulationCase;
        end
        
    end
    
end