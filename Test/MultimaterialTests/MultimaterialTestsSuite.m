classdef MultimaterialTestsSuite < handle

    methods

        function obj = MultimaterialTestsSuite()
            path = './Test/MultimaterialTests/MultimaterialTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','MultiMat');
            results = suite.run;
            table(results)
        end

    end

end