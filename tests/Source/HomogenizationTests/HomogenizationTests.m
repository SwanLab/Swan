classdef HomogenizationTests < testRunner
    
    properties (Access = protected)
        FieldOfStudy = 'Homogenization'
        tests
    end
    
    methods (Access = public)
        function obj = HomogenizationTests()
            obj@testRunner();
        end
    end
    
    methods (Access = protected)
        function loadTests(obj)
            obj.tests = {...
                %% showing error
%                 'testNumericalConvergenceOfNumberOfLaminates';
%                 'testCommutingHomogPlaneStressWithZeroPoisson';
%                 'testCommutingVoigtHomog';
%                 'testAnisotropicPlaneStressbyEnergyEquivalence';
%                 'testStressInPlaneStress';
%                 'testStressRotationInVoigtNotationIn3D';
%                 'testStressRotationInVoigtNotationInPlaneStress';
%                 'testInverseSymmetricFourthOrderTensor' ;
%                 'testInverseNonSymmetricFourthOrderTensor';
%                 'testInverseOfInverseForStiffTensor';
%                 'testIsotropicFourthOrderTensor'
%                 'testInverseSymmetricFourthOrderTensor';
%                 'testSymmetrizeIsotropicFourthOrderTensor';
%                 'testSymmetryForIAniTensorInVoigt'
%                 'testMakeAnisotorpicTensorPlaneStressSymbolically';
%                 'testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor';
%                 'testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor';
%                 'testComplianceTensorThrougtVoigtComparingEnergy';
%                 'TestTwoRankSequentialLaminate';
%                 'testHorizontalTensorRotatedVsVPH';
                %% not showing error
                'testDiagonalLaminate';
                'testNotCommutingHomogPlaneStress';
                'testSymmetrizeFourthOrderTensor';
                'testHorizontalTensorRotatedVsVHP';
                'testHorizontalTensorRotatedVsRank2';
                'testHorizontalTensorRotatedVsHVP';
                'TestGeneralTwoRankSequentialLaminate';
                'testHorizontalLaminate';
                };
        end
        
    end
    
end

