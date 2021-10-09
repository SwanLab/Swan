classdef NewTestSuite < handle
    
    methods

        function obj = NewTestSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
%             suite = matlab.unittest.TestSuite.fromFile('NewFemTests.m');
            suite = matlab.unittest.TestSuite.fromFile('NewFemTests.m', 'Tag','Nou');
            results = suite.run;
            table(results)
        end

    end
end

