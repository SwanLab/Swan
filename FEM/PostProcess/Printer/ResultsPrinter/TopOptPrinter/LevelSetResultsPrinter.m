classdef LevelSetResultsPrinter < ResultsPrinter
    
    properties (Access = protected)
       simulationStr = 'LevelSet';
    end
    
    properties (Access = private)
        fieldName = 'LevelSet';
        headPrinter = NoGaussHeadPrinter;
    end
    
    methods (Access = public)
        
        function obj = LevelSetResultsPrinter(d)
            obj.init(d);
        end
        
    end
    
    methods (Access = protected)
        
        function printHeader(obj)
            obj.headPrinter.print(obj.fileID);
        end
        
        function printResults(obj)
            ls = obj.fields; 
            dS = obj.createScalarDataBase(ls,obj.fieldName,'OnNodes');   
            ScalarPrinter(dS);
        end
    
    end
    

end