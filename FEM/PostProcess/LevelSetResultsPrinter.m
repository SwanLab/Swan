classdef LevelSetResultsPrinter < ResultsPrinter
    
    properties
        fieldName = 'LevelSet';
    end
    
    
    methods (Access = public)
        
        function obj = LevelSetResultsPrinter()
        end
    end
    
    
    methods (Access = protected)
        
        function printHeader(obj)
        end
        
        function printResults(obj)
            ls = obj.fields; 
            iS = obj.istep;
            ScalarPrinter(obj.fileID,ls,obj.fieldName,iS,'OnNodes');            
        end
    
    end
    

end