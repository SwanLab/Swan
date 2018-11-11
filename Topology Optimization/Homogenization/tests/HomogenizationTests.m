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
                               'testDiagonalLaminate';
                                'TestGeneralTwoRankSequentialLaminate';
                               'testHorizontalLaminate';
                               'TestTwoRankSequentialLaminate';
                               'testMakeAnisotrpicTensorPlaneSressSymbollicaly';
                               'testFourthOrderTensor';
                               'testSymmetrizeIsotropicFourthOrderTensor';
                               'testSymmetrizeFourthOrderTensor';
                               'testEnergyEquivalenceVoigtAndTensorNotationForIsoTensor';
                               'testSymmetryForIAniTensorInVoigt';
                               'testEnergyEquivalenceVoigtAndTensorNotationForIAniTensor';
                               'testAnisotropicPlaneStressbyEnergyEquivalence';
                               'testEnergyEquivalenceVoigtAndTensorNotation'};
        end
    end
    
end

