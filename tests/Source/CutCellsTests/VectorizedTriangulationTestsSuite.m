classdef VectorizedTriangulationTestsSuite < handle
    
    methods

        function obj = VectorizedTriangulationTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
            cami = './tests/Source/CutCellsTests/VectorizedTriangulationTests.m'; % des de /swan/
%             cami = 'VectorizedTriangulationTests.m';
            suite = matlab.unittest.TestSuite.fromFile(cami, 'Tag','VectorizedTriangulation'); % VectorizedTriangulation
            results = suite.run;
            table(results)
        end

    end
end