classdef testUnfittedIntegration < testShowingError...
        & testUnfitted
    properties (Access = protected)
        tol = 6e-2;
    end
    
    properties (Access = protected, Abstract)
        testName
        analyticalValue
    end
    
    properties (Access = private)
        varAdim
    end
    
    methods (Access = protected, Abstract)
        
        computeGeometricalVariable(obj)
        
    end
    
    methods (Access = protected)
        
        function obj = testUnfittedIntegration()
            obj.createTopOpt();
            obj.integrateSurface();
        end
        
        function integrateSurface(obj)
            obj.createMesh();
            
            geomVar = obj.computeGeometricalVariable();
            obj.varAdim = geomVar/obj.analyticalValue;
        end
        
        function computeError(obj)
            obj.error = abs(obj.varAdim - 1);
        end
        
    end
    
end

