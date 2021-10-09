classdef NewUnfittedIntegrationTestsSuite < handle
    
    methods

        function obj = NewUnfittedIntegrationTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
%             suite = matlab.unittest.TestSuite.fromFile('NewUnfittedIntegrationTests.m');
            suite = matlab.unittest.TestSuite.fromFile('NewUnfittedIntegrationTests.m', 'Tag','Nou');
            results = suite.run;
            table(results)
        end

    end
end

