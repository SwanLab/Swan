classdef testNumericalHomogenizerPrinter < testNotShowingError...
                                          & testPrintingComparator
    
    properties (Access = protected)
        fileOutputName = 'testNumericalHomogenizerPrinter';
        iter = 0;
    end
    
    properties (Access = private)
        dataBase
    end
    
    methods (Access = public)
        
        function obj = testNumericalHomogenizerPrinter()
            d.iter     = 0;
            d.file     = obj.fileOutputName;      
            d.postCase = 'NumericalHomogenizer';
            tp = testPrintingPostProcessComputer(d);
            tp.print();             
        end
        
    end
   
    
end