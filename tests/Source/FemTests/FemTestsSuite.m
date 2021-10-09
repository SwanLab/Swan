classdef FemTestsSuite < handle
    
    methods

        function obj = FemTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
%             suite = matlab.unittest.TestSuite.fromFile('FemTests.m');
            suite = matlab.unittest.TestSuite.fromFile('FemTests.m', 'Tag','FEM');
            results = suite.run;
            table(results)
        end

    end
end

