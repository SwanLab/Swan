classdef TopOptTestsSuite < handle
    
    methods

        function obj = TopOptTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
            % Nota II: Els Vigdergauz .mat i .m no coincideixen en
            % resultats!
            cami = './tests/Source/TopOptTests/TopOptTests.m'; % des de /swan/
            suite = matlab.unittest.TestSuite.fromFile(cami, 'Tag','TopOpt');
            results = suite.run;
            table(results)
        end

    end
end