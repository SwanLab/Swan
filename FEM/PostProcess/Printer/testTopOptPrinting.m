classdef testTopOptPrinting < testNotShowingError ...
                                  & testTopOptComputation ...
                                  & testPrintingDescriptor
    
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
            obj.init();
            obj.compareFiles();
        end
        
        function init(obj)
            obj.iter = obj.topOpt.settings.maxiter;
        end
        
        function createSettings(obj)
            obj.createSettings@testTopOptComputation();
            obj.settings.printMode = obj.printMode;            
            obj.settings.case_file = obj.fileOutputName;
        end
        
    end
    
    
end

