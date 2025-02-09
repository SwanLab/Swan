classdef SecondOrderTensor2VoigtConverter < Tensor2VoigtConverter...
                                            & Tensor2Voigt3Drepresentation
    
    properties (Access = protected)
        voigtIndex
    end
    
    methods (Access = protected)
        
        function representTensorInVoigt(obj)
            a = obj.tensor.getValue();
            t = obj.voigtTensor.getValue();
            d = obj.tensor.getDimension();
            converter = obj.indexTransformer;
            for i = 1:d
                for j = 1:d
                    iv = converter.tensor2Voigt(i,j);
                    obj.voigtIndex = iv;
                    vf = obj.computeVoigtFactor();
                    t(iv,1) = vf*a(i,j);
                end
            end
            obj.voigtTensor.setValue(t);
        end
        
        function obtainVoigtTensorSize(obj)
            obj.voigtTensorSize = [obj.dimVoigt,1];
        end

    end
    
    methods (Access = protected, Abstract)
        computeVoigtFactor(obj,index)
    end
end

