classdef Tensor2VoigtConverterFactory < TensorVoigtConverterFactory
    
    properties (Access = private)
        tensor
        t2vConverter
    end
    
    methods (Access = public)
        
        function t2vConverter = create(obj,tensor)
            obj.init(tensor)
            obj.createTensor2VoigtConverter()
            t2vConverter = obj.getTensor2VoigtConverter;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,tensor)
            obj.tensor = tensor;
            obj.input = tensor;
        end
          
        function createTensor2VoigtConverter(obj)
            
            if obj.isComplianceTensor()
                
                if obj.isPlaneStress()
                    conv = ComplianceTensor2VoigtConverterPS(obj.tensor);
                elseif obj.is3D()
                    conv = ComplianceTensor2VoigtConverter(obj.tensor);
                else
                    obj.showError();
                end
                
            elseif obj.isStiffnessTensor()
                
                if obj.isPlaneStress()
                    conv = StiffnessTensor2VoigtConverterPS(obj.tensor);
                elseif obj.is3D()
                    conv = StiffnessTensor2VoigtConverter(obj.tensor);
                else
                    obj.showError();
                end
                
                
            elseif obj.isStrainTensor()
                
                if obj.isPlaneStress()
                    conv = StrainTensor2VoigtConverterPS(obj.tensor);
                elseif obj.is3D()
                    conv = StrainTensor2VoigtConverter(obj.tensor);
                else
                    obj.showError();
                end
                
            elseif obj.isStressTensor()
                if obj.isPlaneStress()
                    conv = StressTensor2VoigtConverterPS(obj.tensor);
                elseif obj.is3D()
                    conv = StressTensor2VoigtConverter(obj.tensor);
                else
                    obj.showError();
                end

            else
                obj.showError();
            end
            obj.t2vConverter = conv;
        end
        
        function t2vc = getTensor2VoigtConverter(obj)
            t2vc = obj.t2vConverter;
        end
        
        
    end
    
  
end
