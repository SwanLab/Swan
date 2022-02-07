classdef VectorizedTriangulationTestsSuite < handle
    
    methods

        function obj = VectorizedTriangulationTestsSuite()
            path = './tests/Source/CutCellsTests/VectorizedTriangulationTests.m';
            suite = matlab.unittest.TestSuite.fromFile(path, 'Tag','IsoCoord'); %VectorizedTriangulation, IsoCoord
            results = suite.run;
            table(results)
        end

    end
end