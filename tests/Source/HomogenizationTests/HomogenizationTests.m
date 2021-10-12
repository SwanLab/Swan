classdef HomogenizationTests < handle & matlab.unittest.TestCase

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
            passedTests = {...
                'testDiagonalLaminate';
                'testNotCommutingHomogPlaneStress';
                'testSymmetrizeFourthOrderTensor';
                'testHorizontalTensorRotatedVsVHP';
                'testHorizontalTensorRotatedVsRank2';
                'testHorizontalTensorRotatedVsHVP';
                'TestGeneralTwoRankSequentialLaminate';
                'testHorizontalLaminate';
                }
    end

    methods (Test, TestTags = {'HomogenizationTests', 'ShowingError'})

        function testsError(testCase, errorTests)
            cd ../../../
            test = eval(errorTests);
            err = test.computeError();
            tol = test.tol;
            testCase.verifyLessThanOrEqual(err, tol)
            cd ./tests/Source/HomogenizationTests/
        end

    end

    methods (Test, TestTags = {'HomogenizationTests', 'NotShowingError'})

        function testsPassed(testCase, passedTests)
            cd ../../../
            test = eval(passedTests);
            passed = test.hasPassed();
            verifyTrue(testCase, passed)
            cd ./tests/Source/HomogenizationTests/
        end

    end

end