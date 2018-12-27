classdef DensityResultsPrinter < ResultsPrinter
    
    properties
        fieldName = 'Density';
    end
    
    methods (Access = public)
        
        function obj = DensityResultsPrinter()
        end
    end
    
    
    methods (Access = protected)
        
        function printHeader(obj)
        end
        
        function printResults(obj)
            dens = obj.fields; 
            iS = obj.istep;
            ScalarPrinter(obj.fileID,dens, obj.fieldName,iS,'OnNodes');            
        end
    
    end
    
    
end