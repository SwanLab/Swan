classdef DiffReactTestsSuite < handle

    methods

        function obj = DiffReactTestsSuite()
            path = './Test/DiffReactTests/DiffReactTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','DiffReact');
            results = suite.run;
            table(results)
        end

    end
end