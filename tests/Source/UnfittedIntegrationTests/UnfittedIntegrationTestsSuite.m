classdef UnfittedIntegrationTestsSuite < handle
    
    methods

        function obj = UnfittedIntegrationTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
%             suite = matlab.unittest.TestSuite.fromFile('UnfittedIntegrationTests.m');
            cami = './tests/Source/UnfittedIntegrationTests/UnfittedIntegrationTests.m'; % des de /swan/
%             cami = 'NewUnfittedIntegrationTests.m';
            suite = matlab.unittest.TestSuite.fromFile(cami, 'Tag','UnfittedIntegration');
            results = suite.run;
            table(results)
        end

    end
end

