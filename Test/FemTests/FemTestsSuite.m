classdef FemTestsSuite < handle

    methods

        function obj = FemTestsSuite()
            path = './Test/FemTests/FemTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','FEM');
            results = suite.run;
            table(results)
        end

    end
end