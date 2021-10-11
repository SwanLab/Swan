classdef NewVectorizedTriangulationTestsSuite < handle
    
    methods

        function obj = NewVectorizedTriangulationTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
%             suite = matlab.unittest.TestSuite.fromFile('NewVectorizedTriangulationTests.m');
            cami = './tests/Source/CutCellsTests/NewVectorizedTriangulationTests.m'; % des de /swan/
%             cami = 'NewUnfittedIntegrationTests.m';
            suite = matlab.unittest.TestSuite.fromFile(cami, 'Tag','Nou'); % VectorizedTriangulation
            results = suite.run;
            table(results)
        end

    end
end

