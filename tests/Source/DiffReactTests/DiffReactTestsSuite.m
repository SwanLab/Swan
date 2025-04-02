classdef DiffReactTestsSuite < handle

    methods

        function obj = DiffReactTestsSuite()
            path = './tests/Source/DiffReactTests/DiffReactTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','DiffReact');
            results = suite.run;
            table(results)
        end

    end
end