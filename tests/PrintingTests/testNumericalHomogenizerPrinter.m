classdef testNumericalHomogenizerPrinter < testNotShowingError...
                         & testPrintingDescriptor
    
    properties (Access = protected)
        fileOutputName = 'testNumericalHomogenizerPrinter';
        iter
    end
    
    properties (Access = private)
       dataBase
    end
    
    methods (Access = public)
        
        function obj = testNumericalHomogenizerPrinter()
            obj.iter = 0;
            obj.computeNumericalHomogenizer()
            obj.compareFiles()
        end
        
    end
    
    methods (Access = private)
               
        function computeNumericalHomogenizer(obj)
           dir = Vector3D;
           dir.setValue([1 0 0]);
           lv = 3;
           f = obj.fileOutputName;
           print = true;
           i = obj.iter;
           NumericalFiberHomogenizer(dir,lv,f,print,i);
        end        

    end
    
end