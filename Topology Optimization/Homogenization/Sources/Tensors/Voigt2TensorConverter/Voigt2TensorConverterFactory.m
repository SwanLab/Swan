classdef Voigt2TensorConverterFactory < handle
    
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
        end
        
        function createVoigt2TensorConverter(obj)
            tensorValue = obj.voigtTensor.tensor;
            
            if obj.isFourthOrder()
                
                if obj.isPlaneStress()
                    obj.v2tConverter = PlaneStressFourthOrderVoigt2TensorConverter(tensorValue);
                elseif obj.is3D()
                    obj.v2tConverter = FourthOrderVoigt2TensorConverter(tensorValue);
                else
                    obj.showError();
                end
                
            elseif obj.isStrainTensor()
                
                if obj.isPlaneStress()
                    obj.v2tConverter = StrainVoigt2TensorConverterPS(tensorValue);
                elseif obj.is3D()
                    obj.v2tConverter = StrainVoigt2TensorConverter(tensorValue);
                else
                    obj.showError();
                end
                
            elseif obj.isStressTensor()
                if obj.isPlaneStress()
                    obj.v2tConverter = StressVoigt2TensorConverterPS(tensorValue);
                elseif obj.is3D()
                    obj.v2tConverter = StressVoigt2TensorConverter(tensorValue);
                else
                    obj.showError();
                end
                
            else
                obj.showError();
            end
        end
        
        
        function itIs = isStrainTensor(obj)
            itIs = isa(obj.voigtTensor,'StrainVoigtTensor');
        end
        
        function itIs = isStressTensor(obj)
            itIs = isa(obj.voigtTensor,'StressVoigtTensor');
        end
        
        function itIs = isFourthOrder(obj)
            isVoigt = isa(obj.voigtTensor,'VoigtTensor');
            is3D    = obj.is3D; 
            itIs    = isVoigt || is3D;
        end
        
        function itIs = is3D(obj)            
            itIs = size(obj.voigtTensor.getValue(),1) == 6;
        end
        
        function itIs = isPlaneStress(obj)
            dim = size(obj.voigtTensor.tensor,1);
            itIs = dim == 2;
        end
                
        function t2vc = getVoigt2TensorConverter(obj)
            t2vc = obj.v2tConverter;
        end
        
        
    end
    
    methods (Access = private, Static)
        
        function showError()
            error('Not admitted object to make it tensor')
        end
    end
    
    
    
end

