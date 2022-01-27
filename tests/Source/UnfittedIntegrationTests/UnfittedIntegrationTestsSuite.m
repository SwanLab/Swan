classdef UnfittedIntegrationTestsSuite < handle

    methods

        function obj = UnfittedIntegrationTestsSuite()
            path = './tests/Source/UnfittedIntegrationTests/UnfittedIntegrationTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','UnfittedIntegration');
            results = suite.run;
            table(results)
        end

    end

end