classdef FemTestsSuite < handle

    methods

        function obj = FemTestsSuite()
            path = './tests/Source/FemTests/FemTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','Micro');
            results = suite.run;
            table(results)
        end

    end
end