classdef TestNetwork < matlab.unittest.TestCase
    
        methods(Test)
        function assemblyTest(testCase)
            load("Comprovant.mat","network")
            hiddenlayers    = [1,2];

            data      = Data('../Datasets/Iris.csv',30,1);
            structure = [data.nFeatures,hiddenlayers,data.nLabels];
            actual   = Network(data,structure);
            expected  = network;

            testCase.verifyEqual(actual,expected)
        end
    end
end