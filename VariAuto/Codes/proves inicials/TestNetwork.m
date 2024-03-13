classdef TestNetwork < matlab.unittest.TestCase
    
        methods(Test)
        function assemblyTest(testCase)
            clc;
            load("Comprovant.mat","network")
            hiddenlayers    = [1,2];

            data      =  Data('../Datasets/Iris.csv',30,1);
            structure =  [data.nFeatures,hiddenlayers,data.nLabels];
            ntw       =  Network(data,structure);
            actual    =  optimizationProblem(data,ntw);
            expected  =  network;

            testCase.verifyEqual(actual.thetavec(2),expected.thetavec(2))
        end
    end
end