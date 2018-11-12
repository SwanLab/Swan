classdef Tensor2VoigtConverterFactory < handle
    
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
        end
          
        function createTensor2VoigtConverter(obj)
            tensorValue = obj.tensor.getValue();
            
            if obj.isFourthOrder()
               
                if obj.isPlaneStress()
                    obj.t2vConverter = PlaneStressFourthOrderTensor2VoigtConverter(tensorValue);
                elseif obj.is3D()
                    obj.t2vConverter = FourthOrderTensor2VoigtConverter(obj.tensor);
                else
                    obj.showError();
                end
                
            elseif obj.isStrainTensor()
               
                if obj.isPlaneStress()
                    obj.t2vConverter = StrainTensor2VoigtConverterPS(obj.tensor);
                elseif obj.is3D()
                    obj.t2vConverter = StrainTensor2VoigtConverter(obj.tensor);
                else
                    obj.showError();
                end
                
            elseif obj.isStressTensor()
                if obj.isPlaneStress()
                    obj.t2vConverter = StressTensor2VoigtConverterPS(obj.tensor);
                elseif obj.is3D()
                    obj.t2vConverter = StressTensor2VoigtConverter(obj.tensor);
                else
                    obj.showError();
                end

            else
                obj.showError();
            end
        end
        
        function itIs = isStrainTensor(obj)
            itIs = strcmp(obj.tensor.getFieldName(),'strain');
        end
        
        function itIs = isStressTensor(obj)
            itIs = strcmp(obj.tensor.getFieldName(),'stress');
        end
               
        function itIs = isFourthOrder(obj)
            itIs = strcmp(obj.tensor.getOrder,'fourth');
        end
        
        function itIs = is3D(obj)
            itIs = strcmp(obj.tensor.getElasticityCase,'3D');
        end
        
        function itIs = isPlaneStress(obj)            
            itIs = strcmp(obj.tensor.getElasticityCase,'planeStress');
        end
        
        
        function t2vc = getTensor2VoigtConverter(obj)
            t2vc = obj.t2vConverter;
        end
        
        
    end
    
    methods (Access = private, Static)
        
        function showError()
            error('Not admitted object to make it Voigt')
        end
    end
end
