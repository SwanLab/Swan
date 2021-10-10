classdef NewUnfittedIntegrationTestsSuite < handle
    
    methods

        function obj = NewUnfittedIntegrationTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
%             suite = matlab.unittest.TestSuite.fromFile('NewUnfittedIntegrationTests.m');
            cami = './tests/Source/UnfittedIntegrationTests/NewUnfittedIntegrationTests.m'; % des de /swan/
%             cami = 'NewUnfittedIntegrationTests.m';
            suite = matlab.unittest.TestSuite.fromFile(cami, 'Tag','Nou');
            results = suite.run;
            table(results)
        end

    end
end

