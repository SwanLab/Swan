classdef ZZZTopOptTestsSuite < handle
    
    methods

        function obj = ZZZTopOptTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
%             suite = matlab.unittest.TestSuite.fromFile('FemTests.m');
            cami = './tests/Source/TopOptTests/ZZZTopOptTests.m'; % des de /swan/
%             cami = 'ZZZTopOptTests.m';
            suite = matlab.unittest.TestSuite.fromFile(cami, 'Tag','Nou');
            results = suite.run;
            table(results)
        end

    end
end

