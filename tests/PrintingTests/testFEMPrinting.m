classdef testFEMPrinting < testNotShowingError...
                         & testFemComputation ...
                         & testPrintingDescriptor
    
    properties (Access = protected)
        testName = 'test2d_quad';  
        fileOutputName = 'testFemPrinting';
        iter
    end
    
    properties (Access = private)
       dataBase
    end
    
    methods (Access = public)
        
        function obj = testFEMPrinting()
            obj.iter = obj.fem.getIter();
            obj.print()
            obj.compareFiles()
        end
        
    end
    
    methods (Access = private)
               
        function print(obj)            
            obj.fem.print(obj.fileOutputName);
        end

    end
    
end

