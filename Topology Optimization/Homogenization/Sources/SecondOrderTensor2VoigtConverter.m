classdef SecondOrderTensor2VoigtConverter < Tensor2VoigtConverterFor3DTensors
    
    properties (Access = protected)
        voigtIndex
    end
    
    methods (Access = protected)
        
        function representTensorInVoigt(obj,a)
            t = obj.voigtTensor;
            d = obj.dim;
            converter = obj.indexTransformer; 
            for i = 1:d
                for j = 1:d
                    iv = converter.tensor2Voigt(i,j);
                    obj.voigtIndex = iv;
                    vf = obj.computeVoigtFactor();
                    t(iv) = vf*a(i,j);
                end
            end
            obj.voigtTensor = t;
        end
        
        function obtainVoigtTensorSize(obj)
            obj.voigtTensorSize = [obj.dimVoigt,1];
        end

    end
    methods (Access = protected, Abstract)
        computeVoigtFactor(obj,index)
    end
end

