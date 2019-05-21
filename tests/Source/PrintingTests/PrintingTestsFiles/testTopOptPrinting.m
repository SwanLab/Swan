classdef testTopOptPrinting < testNotShowingError ...
                            & testPrintingComparator 
                                  
    
    properties (Access = protected)
      iter
      fields
      dataBase
    end
       
    properties (Access = protected, Abstract)
        fileOutputName
        printMode
    end
    
    methods (Access = protected)
        
        function obj = testTopOptPrinting()
            obj.iter = 0;
            d.iter     = obj.iter;
            d.file     = obj.fileOutputName;      
            d.postCase = 'TopOptProblem';
            tp = testPrintingPostProcessComputer(d);
            tp.print();
        end
         
    end    
    
end

