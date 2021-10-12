classdef ImageProcessingTestsSuite < handle
    
    methods

        function obj = ImageProcessingTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
            suite = matlab.unittest.TestSuite.fromFile('ImageProcessingTests.m', 'Tag','ImageProcessing');
            results = suite.run;
            table(results)
        end

    end
end