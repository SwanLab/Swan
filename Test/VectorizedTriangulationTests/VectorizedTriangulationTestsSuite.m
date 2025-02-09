classdef VectorizedTriangulationTestsSuite < handle
    
    methods

        function obj = VectorizedTriangulationTestsSuite()
            path = './Test/VectorizedTriangulationTests/VectorizedTriangulationTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','VectorizedTriangulation');
            results = suite.run;
            table(results)
        end

    end
end