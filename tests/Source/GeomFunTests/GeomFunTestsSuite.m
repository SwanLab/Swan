classdef GeomFunTestsSuite < handle

    methods

        function obj = GeomFunTestsSuite()
            path = './tests/Source/GeomFunTests/GeomFunTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','GeomFun');
            results = suite.run;
            table(results)
        end

    end

end