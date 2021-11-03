classdef TestPlotting < testUnfitted
    
    properties (Access = protected)
        testName;
        meshType;
        meshIncludeBoxContour;
    end

    properties (Access = private)
        varAdim
    end
    
    methods (Access = public)

        
        function obj = TestPlotting(cParams)
            obj.init(cParams);
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

        function passed = computePassed(obj)
            d = load(obj.testName);
            itIs = isequaln(obj.unfittedMesh,d.unfittedMesh);
            passed = itIs;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.testName        = cParams.testName;
            obj.meshType        = cParams.meshType;
            obj.meshIncludeBoxContour = cParams.meshIncludeBoxContour;
        end

        function angle = getViewAngle(obj)
            if isprop(obj,'viewAngle')
                angle = obj.viewAngle;
            else
                angle = [0 0 1];
            end
        end
        
    end

    %% Heredat de testUnfitted, s'ha de declarar buit
    methods (Access = protected)
        function printTestNotPassed()
        end
        function printTestPassed()
        end
        function hasPassed()
        end
    end
end