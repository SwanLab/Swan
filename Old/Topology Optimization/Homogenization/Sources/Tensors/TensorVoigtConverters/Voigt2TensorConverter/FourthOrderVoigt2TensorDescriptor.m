classdef FourthOrderVoigt2TensorDescriptor < handle
    
    properties (Access = protected, Abstract) 
        tensor
        indexTransformer
        voigtTensor
    end
    
    
    methods (Access = protected)
        
       function  representVoigtInTensor(obj)
            a = obj.voigtTensor.getValue();
            c = obj.tensor.getValue();
            converter = obj.indexTransformer;
            d  = obj.voigtTensor.getVoigtDimension();
            for iv = 1:d
                for jv = 1:d
                    [i,j] = converter.voigt2tensor(iv);
                    [k,l] = converter.voigt2tensor(jv);
                    factor = obj.getVoigtFactor(iv,jv);
                    aij = 1/factor*a(iv,jv);
                    
                    c(i,j,k,l) = aij;
                    c(j,i,k,l) = aij;
                    c(i,j,l,k) = aij;
                    c(j,i,l,k) = aij;
                    
                    c(k,l,i,j) = aij;
                    c(l,k,i,j) = aij;
                    c(k,l,j,i) = aij;
                    c(l,k,j,i) = aij;
                end
            end
            obj.tensor.setValue(c);
        end
        
    end
    
    methods (Access = protected, Abstract, Static)
        getVoigtFactor()
    end
    
end