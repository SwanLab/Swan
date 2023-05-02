classdef TestComputerIPM < handle

    properties (Access = public)
    end
    properties (Access = private)  
        tolerance
        maxDifference
        desiredTest
        actualData
        loadedData 
    end

    methods (Access = public)
        function obj = TestComputerIPM(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeMaxDifference();
            obj.checkStatus();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.actualData = cParams.actualData;
            obj.loadedData = cParams.loadedData;
            obj.desiredTest = cParams.desiredTest;
            obj.tolerance = 1e-5;
        end

        function computeMaxDifference(obj)
            obj.maxDifference = abs(max(max(obj.loadedData - obj.actualData.value')));
        end

        function checkStatus(obj)
            if obj.maxDifference > obj.tolerance
                error('The following test has failed: %s Test',obj.desiredTest);
            else
                fprintf(obj.desiredTest); fprintf(' test passed \n');
            end
        end
    end
end