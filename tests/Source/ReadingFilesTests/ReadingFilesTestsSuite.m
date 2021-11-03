classdef ReadingFilesTestsSuite < handle
    
    methods

        function obj = ReadingFilesTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
            suite = matlab.unittest.TestSuite.fromFile('ReadingFilesTests.m', 'Tag','ReadingFiles');
            results = suite.run;
            table(results)
        end

    end
end