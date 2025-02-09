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
            elementdiffreact = {'testHorizontalTensorRotatedVsHVP'};
    end

%     methods (Test, TestTags = {'HomogenizationTests', 'Elementdiffreact'})
% 
%         function testsElementdiffreact(testCase, elementdiffreact)
%             testCase.fixFolder();
%             test = eval(elementdiffreact);
%             passed = test.hasPassed();
%             verifyTrue(testCase, passed)
%         end
% 
%     end

    methods (Test, TestTags = {'HomogenizationTests', 'ShowingError'})

        function testsError(testCase, errorTests)
            testCase.fixFolder();
            test = eval(errorTests);
            err = test.computeError();
            tol = test.tol;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Test, TestTags = {'HomogenizationTests', 'NotShowingError'})

        function testsPassed(testCase, passedTests)
            testCase.fixFolder();
            test = eval(passedTests);
            passed = test.hasPassed();
            verifyTrue(testCase, passed)
        end

    end

    methods (Access = private)
        
        function fixFolder(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            changeToFolder = '../../';
            testCase.applyFixture(CurrentFolderFixture(changeToFolder));
        end
    end

end