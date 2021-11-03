classdef PlottingTestsSuite < handle
    
    methods

        function obj = PlottingTestsSuite()
            % Nota: queda corregir que els tests puguin correr fora
            % d'aquest directori
            close all;
            cami = './tests/Source/PlottingTests/PlottingTests.m'; % des de /swan/
            suite = matlab.unittest.TestSuite.fromFile(cami, 'Tag','PlottingTests');
            results = suite.run;
            table(results)
        end

    end
end

