classdef ConvectionDiffusionTestsSuite < handle

    methods

        function obj = ConvectionDiffusionTestsSuite()
            path = './HypersonicAeroRawFiles/SteadyConvectionDiffusion/ConvectionDiffusionTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','ConvectionDiffusion');
            results = suite.run;
            table(results)
        end

    end

end