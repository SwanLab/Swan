classdef ZZPlottingTestsSuite < handle
    
    methods

        function obj = ZZPlottingTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
            cami = './tests/Source/PlottingTests/ZZPlottingTests.m'; % des de /swan/
            suite = matlab.unittest.TestSuite.fromFile(cami, 'Tag','PlottingTests');
            results = suite.run;
            table(results)
        end

    end
end

