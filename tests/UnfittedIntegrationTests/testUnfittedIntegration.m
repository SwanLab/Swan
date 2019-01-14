classdef testUnfittedIntegration < testShowingError...
                                   & testUnfitted
    properties (Access = protected)
        tol = 6e-2;
    end
    
    properties (Access = protected, Abstract)
        testName
        analiticalArea
    end
    
    properties (Access = private)
        areaAdim
    end
    
    methods (Access = protected)
        function obj = testUnfittedIntegration()
            obj.createTopOpt()
            obj.integrateSurface()
        end
        
        function integrateSurface(obj)
            obj.createMesh();
            
            area = obj.computeGeometricalVariable();
            obj.areaAdim = area/obj.analiticalArea;
        end
        
        function M = computeGeometricalVariable(obj)
            M = obj.mesh.computeMass();
        end
        
        function computeError(obj)
            obj.error = abs(obj.areaAdim - 1);
        end
    end
end

