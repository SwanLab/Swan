classdef ZZTestPlotting < testUnfitted
    
    properties (Access = protected)
        testName;
        meshType;
        meshIncludeBoxContour;
    end

    properties (Access = private)
        varAdim
    end
    
    methods (Access = public)

        
        function obj = ZZTestPlotting(cParams)
            obj.init(cParams);
            obj.createTopOpt(); % peta aqui
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

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.testName        = cParams.testName;
            obj.meshType        = cParams.meshType;
            obj.meshIncludeBoxContour = cParams.meshIncludeBoxContour;
        end
        
    end

    methods (Access = public)

        function passed = computeError(obj)
            d = load(obj.testName);
            itIs = isequaln(obj.unfittedMesh,d.unfittedMesh);
            passed = itIs;
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