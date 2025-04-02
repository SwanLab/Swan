classdef FourthOrderTensor2VoigtDescriptor < handle
    
    properties (Access = protected, Abstract)
        tensor
        indexTransformer
        voigtTensor
    end
    
    
    methods (Access = protected)
        
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
                            factor = obj.getVoigtFactor(iv,jv);
                            c(iv,jv) = factor*a(i,j,k,l);
                        end
                    end
                end
            end
            obj.voigtTensor.setValue(c);
        end
        
    end
    
    methods (Access = protected, Abstract, Static)
        getVoigtFactor()
    end


end