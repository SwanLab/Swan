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
            tensorValue = obj.tensor.tensor;
            
            if obj.isFourthOrder()
               
                if obj.isPlaneStress()
                    obj.t2vConverter = PlaneStressFourthOrderTensor2VoigtConverter(tensorValue);
                elseif obj.is3D()
                    obj.t2vConverter = FourthOrderTensor2VoigtConverter(tensorValue);
                else
                    obj.showError();
                end
                
            elseif obj.isStrainTensor()
               
                if obj.isPlaneStress()
                    obj.t2vConverter = StrainTensor2VoigtConverterPS(tensorValue);
                elseif obj.is3D()
                    obj.t2vConverter = StrainTensor2VoigtConverter(tensorValue);
                else
                    obj.showError();
                end
                
            elseif obj.isStressTensor()
                if obj.isPlaneStress()
                    obj.t2vConverter = StressTensor2VoigtConverterPS(tensorValue);
                elseif obj.is3D()
                    obj.t2vConverter = StressTensor2VoigtConverter(tensorValue);
                else
                    obj.showError();
                end

            else
                obj.showError();
            end
        end
        
        function itIs = isStrainTensor(obj)
            itIs = isa(obj.tensor,'StrainTensor');
        end
        
        function itIs = isStressTensor(obj)
            itIs = isa(obj.tensor,'StressTensor');
        end
               
        function itIs = isFourthOrder(obj)
            first  = isa(obj.tensor,'FourthOrderTensor');
            second = isa(obj.tensor,'IsotropicConstitutiveTensor3D');
            itIs = first || second;
        end
        
        function itIs = is3D(obj)
            dim = size(obj.tensor.tensor,1);
            itIs = dim == 3;
        end
        
        function itIs = isPlaneStress(obj)
            dim = size(obj.tensor.tensor,1);
            itIs = dim == 2;
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
