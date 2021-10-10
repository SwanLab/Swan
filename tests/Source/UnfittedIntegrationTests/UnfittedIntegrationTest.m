classdef UnfittedIntegrationTest < testUnfitted
    
    properties (Access = protected)
        testName;
        analyticalValue;
        meshType;
        meshIncludeBoxContour;
    end

    properties (Access = private)
        varAdim
    end

    %% Heredat de testUnfittedPerimeterIntegration
    
    methods (Access = public)

        
        function obj = UnfittedIntegrationTest(cParams)
            obj.init(cParams);
            obj.createTopOpt();
            obj.integrateSurface();
        end
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.testName = cParams.testName;
            obj.analyticalValue = cParams.analyticalValue;
            obj.meshType = cParams.meshType;
            obj.meshIncludeBoxContour = cParams.meshIncludeBoxContour;
        end
        
        function integrateSurface(obj)
            obj.createMesh();
            geomVar = obj.computeGeometricalVariable();
            obj.varAdim = geomVar/obj.analyticalValue;
        end
        
        function totalIntegral = computeGeometricalVariable(obj)
            switch obj.meshType  
                case 'INTERIOR'
                    totalIntegral = obj.unfittedMesh.computeMass();
                case 'BOUNDARY'
                    totalIntegral = obj.unfittedMesh.computePerimeter();
            end

        end
        
    end

    %% Heredat de testUnfittedIntegration
    methods (Access = public)
        
        function error = computeError(obj)
            error = abs(obj.varAdim - 1);
        end
        
    end

    methods (Access = protected)
        function printTestNotPassed()
        end
        function printTestPassed()
        end
        function hasPassed()
        end
    end
end

