classdef FourthOrderTensor2VoigtConverter < Tensor2VoigtConverterFor3DTensors
    

    methods (Access = public)
        
        function obj = FourthOrderTensor2VoigtConverter(Tensor)
            obj.computeConversion(Tensor)
        end        
    end
    
    methods (Access = protected)
        
        
        function selectVoigtTensorClass(obj)
            if obj.tensor.getElasticityCase == '3D'
                obj.voigtTensor = SymmetricFourthOrder3DVoigtTensor();
            elseif obj.tensor.getElasticityCase == 'PlaneStress'
                obj.voigtTensor = SymmetricFourthOrderPlaneStressVoigtTensor();
            end
        end
        
       function  representTensorInVoigt(obj)
            a = obj.tensor.getValue();
            c = obj.voigtTensor.getValue();
            converter = obj.indexTransformer; 
            d  = obj.tensor.getDimension();
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l = 1:d
                            iv = converter.tensor2Voigt(i,j);
                            jv = converter.tensor2Voigt(k,l);
                            c(iv,jv) = a(i,j,k,l);
                        end
                    end
                end
            end
            obj.voigtTensor.setValue(c);
       end
        
       
    end
    
end

