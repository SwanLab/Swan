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
                'testNumericalConvergenceOfNumberOfLaminates'
                'testDiagonalLaminate';
                'TestGeneralTwoRankSequentialLaminate';
                'testHorizontalLaminate';
                'TestTwoRankSequentialLaminate';
                'testMakeAnisotorpicTensorPlaneStressSymbolically';
                'testIsotropicFourthOrderTensor';
                'testSymmetrizeIsotropicFourthOrderTensor';
                'testSymmetryForIAniTensorInVoigt'
                'testStressInPlaneStress'
                'testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor';
                'testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor';
                'testAnisotropicPlaneStressbyEnergyEquivalence'
                'testSymmetrizeFourthOrderTensor'};

        end
    end
    
end

