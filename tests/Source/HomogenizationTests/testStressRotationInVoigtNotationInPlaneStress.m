classdef testStressRotationInVoigtNotationInPlaneStress < testStressRotationInVoigtNotation

    methods (Access = public)
        
        function obj = testStressRotationInVoigtNotationInPlaneStress()
            obj.compute()
        end
        
    end
    
    methods (Access = protected)

        function createDirection(obj)
            dir = [ 0 0 1];
            obj.direction = Vector3D;
            obj.direction.setValue(dir);
        end

        function createStress(obj)
            obj.createStress@testStressRotationInVoigtNotation()
            sv = obj.stressVoigt;
            obj.stressVoigt = PlaneStressTransformer.transform(sv);
        end

        function createRotatedStress(obj)
            obj.createRotatedStress@testStressRotationInVoigtNotation()
            sR = obj.rotatedStress;
            obj.rotatedStress = PlaneStressTransformer.transform(sR);
        end
 
    end

end