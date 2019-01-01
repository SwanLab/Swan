classdef testTopOptLevelSetPrinting < testTopOptPrinting
    
    properties (Access = protected)
        testName = 'test_gripping';  
        fileOutputName = 'testTopOptLevelSetPrinting';
        postProcessor = 'NodalLevelSet';
    end
    
        methods (Access = protected)
        
        function computeFields(obj)
            obj.fields = obj.topOpt.x;
        end
        
    end
    
end

