classdef AAHomogenizationTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
            errorTests = {...
                'testNumericalConvergenceOfNumberOfLaminates';
                'testCommutingHomogPlaneStressWithZeroPoisson';
                'testCommutingVoigtHomog';
                'testAnisotropicPlaneStressbyEnergyEquivalence';
                'testStressInPlaneStress';
                'testStressRotationInVoigtNotationIn3D';
                'testStressRotationInVoigtNotationInPlaneStress';
                'testInverseSymmetricFourthOrderTensor' ;
                'testInverseNonSymmetricFourthOrderTensor';
                'testInverseOfInverseForStiffTensor';
                'testIsotropicFourthOrderTensor'
                'testInverseSymmetricFourthOrderTensor';
                'testSymmetrizeIsotropicFourthOrderTensor';
                'testSymmetryForIAniTensorInVoigt'
                'testMakeAnisotorpicTensorPlaneStressSymbolically';
                'testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor';
                'testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor';
                'testComplianceTensorThrougtVoigtComparingEnergy';
                'TestTwoRankSequentialLaminate';
                'testHorizontalTensorRotatedVsVPH';
                }
    end

    methods (Test, TestTags = {'HomogenizationTests', 'ShowingError', 'Nou'})

        function testsError(testCase, errorTests)
            cd ../../../
            test = eval(errorTests);
            err = test.computeError();
            tol = test.tol;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/HomogenizationTests/
        end

    end

end