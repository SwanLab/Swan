classdef HomogenizationTests < testRunner
    
    
    properties (Access = protected)
        FieldOfStudy = 'TensorAndVoigtNotation'
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
               %'test2DSeqLaminateInVoigtWithFormuleOfAllaireWebPAgeExercise'
               'testStressRotationInVoigtNotationIn3D';
                'testStressRotationInVoigtNotationInPlaneStress';               
                'testInverseSymmetricFourthOrderTensor' ;
                'testHorizontalLaminate';
                'testAnisotropicPlaneStressbyEnergyEquivalence';
                'testCommutingVoigtHomogPlaneStress';
                'testInverseOfInverseForStiffTensor';
                'testComplianceTensorThrougtVoigtComparingEnergy';
                'testIsotropicFourthOrderTensor'
                'testInverseSymmetricFourthOrderTensor';
                'testSymmetrizeIsotropicFourthOrderTensor';
                'testSymmetryForIAniTensorInVoigt'
                'testInverseFourthOrderTensor';
                'testMakeAnisotorpicTensorPlaneStressSymbolically';
                'testStressInPlaneStress'
                'testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor';
                'testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor';
                'testSymmetrizeFourthOrderTensor';
                'TestTwoRankSequentialLaminate';
                'testNumericalConvergenceOfNumberOfLaminates';
                'testDiagonalLaminate';
                'TestGeneralTwoRankSequentialLaminate';
                'testHorizontalTensorRotatedVsVHP';
                'testHorizontalTensorRotatedVsVPH';                
                'testHorizontalTensorRotatedVsRank2';
                'testHorizontalTensorRotatedVsHVP';                
                };

        end
    end
    
end

