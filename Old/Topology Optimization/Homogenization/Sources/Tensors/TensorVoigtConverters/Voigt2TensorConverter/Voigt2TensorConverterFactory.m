classdef Voigt2TensorConverterFactory < TensorVoigtConverterFactory
    
    properties (Access = private)
        voigtTensor
        v2tConverter
    end
    
    methods (Access = public)
        
        function t2vConverter = create(obj,tensor)
            obj.init(tensor)
            obj.createVoigt2TensorConverter()
            t2vConverter = obj.getVoigt2TensorConverter();
        end
    end
    
    methods (Access = private)
        
        function init(obj,tensor)
            obj.voigtTensor = tensor;
            obj.input = tensor;
        end
        
        function createVoigt2TensorConverter(obj)
            tensorValue = obj.voigtTensor.getValue();
            
            if obj.isComplianceTensor()
                
                if obj.isPlaneStress()
                    conv = ComplianceVoigt2TensorConverterPS(obj.voigtTensor);
                elseif obj.is3D()
                    conv = ComplianceVoigt2TensorConverter(obj.voigtTensor);
                else
                    obj.showError();
                end
                
            elseif obj.isStiffnessTensor()
                
                if obj.isPlaneStress()
                    conv = StiffnessVoigt2TensorConverterPS(obj.voigtTensor);
                elseif obj.is3D()
                    conv = StiffnessVoigt2TensorConverter(obj.voigtTensor);
                else
                    obj.showError();
                end
                
                
            elseif obj.isStrainTensor()
                
                if obj.isPlaneStress()
                    conv = StrainVoigt2TensorConverterPS(tensorValue);
                elseif obj.is3D()
                    conv = StrainVoigt2TensorConverter(tensorValue);
                else
                    obj.showError();
                end
                
            elseif obj.isStressTensor()
                
                if obj.isPlaneStress()
                    conv = StressVoigt2TensorConverterPS(tensorValue);
                elseif obj.is3D()
                    conv = StressVoigt2TensorConverter(tensorValue);
                else
                    obj.showError();
                end
                
            else
                obj.showError();
            end
            
            obj.v2tConverter = conv;
        end   
        
        function t2vc = getVoigt2TensorConverter(obj)
            t2vc = obj.v2tConverter;
        end
        
    end
    
end

