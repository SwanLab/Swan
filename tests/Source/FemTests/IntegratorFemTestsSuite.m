classdef IntegratorFemTestsSuite < handle
    
    methods

        function obj = IntegratorFemTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
            suite = matlab.unittest.TestSuite.fromFile('IntegratorFemTests.m');
%             suite = matlab.unittest.TestSuite.fromFile('IntegratorFemTests.m', 'Tag','Triangle');
            results = suite.run;
            table(results)
        end

    end
end

