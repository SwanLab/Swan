classdef testNumericalHomogenizerPrinter < testNotShowingError...
        & testPrintingDescriptor
    
    properties (Access = protected)
        fileOutputName = 'testNumericalHomogenizerPrinter';
        microFile = 'test2d_micro';
        iter = 0;
    end
    
    properties (Access = private)
        dataBase
    end
    
    methods (Access = public)
        
        function obj = testNumericalHomogenizerPrinter()
            d = obj.createNumericalHomogenizerDataBase();
            NumericalHomogenizer(d);
        end
        
    end
    
    methods (Access = private)
        
        function d = createNumericalHomogenizerDataBase(obj)
            nDB = NumericalHomogenizerDataBase(obj.microFile);
            d = nDB.dataBase;
            d.print = true;
            d.outFileName = obj.fileOutputName;            
        end
    end
    
        
    
end