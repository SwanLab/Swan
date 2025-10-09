classdef AlgebraicOperationsTestsSuite < handle

    methods

        function obj = AlgebraicOperationsTestsSuite()
            path = './Test/AlgebraicOperationsTests/AlgebraicOperationsTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','Algebra');
            results = suite.run;
            table(results)
        end

    end

end