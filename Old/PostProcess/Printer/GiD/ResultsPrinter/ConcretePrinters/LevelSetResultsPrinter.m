classdef LevelSetResultsPrinter < DesignVariablePrinter
    
    properties (Access = protected)
        simulationStr = 'LevelSet';
        fieldName = 'LevelSet';
    end
    
    methods (Access = public)
        
        function obj = LevelSetResultsPrinter(d)
            obj.init(d);
        end
        
    end
    
end

