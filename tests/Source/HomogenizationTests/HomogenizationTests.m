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
                'testHorizontalLaminate';
                'TestGeneralTwoRankSequentialLaminate';
                'testNotCommutingHomogPlaneStress';
                'testSymmetrizeFourthOrderTensor';
                'testHorizontalTensorRotatedVsVHP';
                'testHorizontalTensorRotatedVsRank2';
                'testHorizontalTensorRotatedVsHVP';
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