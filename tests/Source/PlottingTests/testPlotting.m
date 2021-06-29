classdef testPlotting < testNotShowingError...
        & testUnfitted...
        
    properties (Access = protected, Abstract)
        testName
        meshType
    end
    
    methods (Access = protected)
        
        function obj = testPlotting()
            obj.createTopOpt();
            obj.createMesh();
            obj.plot();
        end
        
        function plot(obj)
            figure();
            if isequal(obj.meshType,'BOUNDARY')
                obj.unfittedMesh.plotBoundary();
                obj.unfittedMesh.plotNormals();
            elseif isequal(obj.meshType,'INTERIOR')
                obj.unfittedMesh.plot();
            end
            view(obj.getViewAngle());
        end
        
        function hasPassed = hasPassed(obj)
            d = load(obj.testName);
            itIs = isequaln(obj.unfittedMesh,d.unfittedMesh);
            hasPassed = itIs;
        end
        
    end
    
    methods (Access = private)
        

        function angle = getViewAngle(obj)
            if isprop(obj,'viewAngle')
                angle = obj.viewAngle;
            else
                angle = [0 0 1];
            end
        end
    end
end

