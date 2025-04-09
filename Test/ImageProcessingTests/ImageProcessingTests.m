classdef ImageProcessingTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        imageTests = { 'test_DenoisingEinstein'; ...
            'test_AcceleratedDenoisingEinstein'}
    end

    methods (Test, TestTags = {'ImageProcessing', 'Fast', 'Gmsh'})

        function testDisplacement(testCase, imageTests)
            run(imageTests);
            s.imageFile            = imageFile;
            s.lipschitzConstant    = lipschitzConstant;
            s.totalVariationWeight = totalVariationWeight;
            s.noiseAmplitud        = noiseAmplitud;
            s.maxIter              = maxIter;
            s.optimizer            = optimizer;
            s.computerType         = 'IMAGE';
            s.testName             = imageTests;
            s.variablesToStore     = {'optimizedImage'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 5e-2;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end
    
end