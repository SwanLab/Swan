classdef HomogenizationTestsSuite < handle
    
    methods

        function obj = HomogenizationTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
            cami = './tests/Source/HomogenizationTests/HomogenizationTests.m'; % des de /swan/
            suite = matlab.unittest.TestSuite.fromFile(cami, 'Tag','HomogenizationTests');
            results = suite.run;
            table(results)
        end

    end
end