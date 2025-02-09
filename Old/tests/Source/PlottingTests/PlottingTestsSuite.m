classdef PlottingTestsSuite < handle
    
    methods

        function obj = PlottingTestsSuite()
            path = './tests/Source/PlottingTests/PlottingTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','PlottingTests');
            results = suite.run;
            table(results)
        end

    end
end