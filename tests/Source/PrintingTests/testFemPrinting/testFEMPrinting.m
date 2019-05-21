classdef testFEMPrinting < testNotShowingError...
                           & testPrintingComparator
                         
    
    properties (Access = protected)
        fileOutputName = 'testFemPrinting';
        iter
    end
    
    properties (Access = private)
       dataBase
    end
    
    methods (Access = public)
        
        function obj = testFEMPrinting() 
            obj.iter   = 0;
            d.iter     = obj.iter;            
            d.file     = obj.fileOutputName;      
            d.postCase = 'Elasticity';
            tp = testPrintingPostProcessComputer(d);
            tp.print();
        end
        
    end
    
end

