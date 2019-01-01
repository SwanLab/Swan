classdef testTopOptDensityPrinting < testTopOptPrinting
    
    properties (Access = protected)
        testName = 'test_cantilever';  
        fileOutputName = 'testTopOptDensityPrinting';
        postProcessor = 'NodalDensity';
    end
    
    
    methods (Access = protected)
        
        function computeFields(obj)
            obj.fields = obj.topOpt.x;
        end
        
    end

    
    
end

