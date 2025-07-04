classdef HyperelasticityTestsSuite < handle

    methods

        function obj = HyperelasticityTestsSuite()
            path = './Test/HyperelasticityTests/HyperelasticityTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','DisplReact');
            results = suite.run;
            table(results)
        end

    end

end