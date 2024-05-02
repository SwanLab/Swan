classdef UnsteadyConvectionTestsSuite < handle

    methods

        function obj = UnsteadyConvectionTestsSuite()
            path = './HypersonicAeroFiles/UnsteadyConvection/Tests/UnsteadyConvectionTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','UnsteadyConvection');
            results = suite.run;
            table(results)
        end

    end

end