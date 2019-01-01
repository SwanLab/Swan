classdef LevelSetResultsPrinter < ResultsPrinter
    
    properties (Access = private)
        fieldName = 'LevelSet';
        simulationCase = 'LevelSet';
        headPrinter = NoGaussHeadPrinter;
    end
    
    
    methods (Access = public)
        
        function obj = LevelSetResultsPrinter()
        end
    end
    
    
    methods (Access = protected)
        
        function printHeader(obj)
            obj.headPrinter.print(obj.fileID);
        end
        
        function printResults(obj)
            ls = obj.fields; 
            iS = obj.istep;
            dS = obj.createScalarDataBase(obj.fileID,ls,obj.fieldName,iS,'OnNodes');   
            ScalarPrinter(dS);
        end
    
    end
    
    methods (Access = private)
        
        function d = createScalarDataBase(obj,fileID, fieldValues,fieldName, istep, fieldPosition)
            d.fileID = fileID;
            d.fieldValues = fieldValues;
            d.fieldName = fieldName;
            d.istep = istep;
            d.fieldPosition = fieldPosition;
            d.simulationCase = obj.simulationCase;
        end
        
    end
    

end