classdef DensityResultsPrinter < DesignVariablePrinter
    
    
    properties (Access = protected)
        simulationStr = 'NodalDensity';
        fieldName   = 'Density';
    end
    
    methods (Access = public)
        
        function obj = DensityResultsPrinter(d)
            obj.init(d);
        end
        
    end
    
end