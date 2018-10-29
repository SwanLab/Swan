classdef testStressRotationInVoigtNotationInPlaneStress < testStressRotationInVoigtNotation
    
    properties (Access = private) 
    end
    
    methods (Access = public)
        
        function obj = testStressRotationInVoigtNotationInPlaneStress()
            obj.compute()
        end
        
    end
    
    methods (Access = protected)
        
        function createDirection(obj)
            obj.Direction = [0 0 1];
        end    
        
        function createStress(obj)
            obj.createStress@testStressRotationInVoigtNotation()
            sv = obj.StressVoigt.getValue();
            tensPS = PlaneStressTransformer.transform(sv);
            obj.StressVoigt.setValue(tensPS)
        end
        
        function createRotatedStress(obj)
            obj.createRotatedStress@testStressRotationInVoigtNotation()
            tens = obj.rotatedStress.getValue();
            tensorPS = PlaneStressTransformer.transform(tens);
            obj.rotatedStress.setValue(tensorPS)
        end
            
    end
    
        
    
end