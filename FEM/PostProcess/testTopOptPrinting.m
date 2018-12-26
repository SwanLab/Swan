classdef testTopOptPrinting < testNotShowingError ...
                                  & testTopOptComputation ...
                                  & testPrintingDescriptor
    
    properties (Access = protected)
      filesHaveChanged
      iter
    end
    
    properties (Access = protected, Abstract)
       postProcessor 
    end
    
    methods (Access = protected)
        
        function obj = testTopOptPrinting()
            obj.init();
            obj.print();
            obj.compareFiles();
        end
        
        function init(obj)
            obj.iter = 0;
        end
        
        function print(obj)
            field     = obj.topOpt.x;
            outName   = obj.fileOutputName;
            mesh      = obj.topOpt.mesh;
            it        = obj.iter;
            postprocess = Postprocess.create(obj.postProcessor);
            postprocess.print(mesh,field,it,outName);
        end
        
    end
    
end

