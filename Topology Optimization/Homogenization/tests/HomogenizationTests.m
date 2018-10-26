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
                'testHorizontalLaminate';
                'testNumericalConvergenceOfNumberOfLaminates';
                'testDiagonalLaminate';
                'TestGeneralTwoRankSequentialLaminate';
                'testHorizontalTensorRotatedVsRank2';
                'testHorizontalTensorRotatedVsVHP';
                'testHorizontalTensorRotatedVsHVP';
                'testHorizontalTensorRotatedVsVPH';
                'testCommutingVoigtHomogPlaneStress';
                'testInverseOfInverseForStiffTensor';
                'testComplianceTensorThrougtVoigtComparingEnergy';
                'testIsotropicFourthOrderTensor'
                'testInverseSymmetricFourthOrderTensor';
                'testStressRotationInVoigtNotationInPlaneStress';
                'testSymmetrizeIsotropicFourthOrderTensor';
                'testSymmetryForIAniTensorInVoigt'
                'testInverseFourthOrderTensor';
                'testStressRotationInVoigtNotationIn3D';
                'testMakeAnisotorpicTensorPlaneStressSymbolically';
                'testStressInPlaneStress'
                'testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor';
                'testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor';
                'testAnisotropicPlaneStressbyEnergyEquivalence'
                'testSymmetrizeFourthOrderTensor';
                'TestTwoRankSequentialLaminate';
};

        end
    end
    
end

