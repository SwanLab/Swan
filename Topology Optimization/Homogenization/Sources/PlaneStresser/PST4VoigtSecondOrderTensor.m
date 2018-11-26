classdef PST4VoigtSecondOrderTensor < PlaneStressTransformer
    
    
    properties (Access = private)
        tensor
    end
    
    methods (Access = public)
        
        function obj = PST4VoigtSecondOrderTensor(tensor)
            obj.tensor = tensor;
            obj.createPlaneStressTensor()
            obj.computePlaneStressTensor()
        end
        
    end
    
    methods (Access = private)
        
        function itIs = isStrain(obj)
            itIs = strcmp(obj.tensor.getFieldName,'strain');
        end
        
        function itIs = isStress(obj)
            itIs = strcmp(obj.tensor.getFieldName,'stress');
        end
        
        function computePlaneStressTensor(obj)
            PSIndex = PlaneStressIndex();
            InPlaneIndex = PSIndex.getInPlaneIndex();
            tens = obj.tensor.getValue;
            tensPS = tens(InPlaneIndex);
            obj.psTensor.setValue(tensPS);
        end
        
    end
    
    methods (Access = protected)
        function createPlaneStressTensor(obj)
            if obj.isStrain()
                obj.psTensor = StrainPlaneStressVoigtTensor();
            elseif obj.isStress()
                obj.psTensor = StressPlaneStressVoigtTensor();
            end
            
        end
        
    end

  
    
end

