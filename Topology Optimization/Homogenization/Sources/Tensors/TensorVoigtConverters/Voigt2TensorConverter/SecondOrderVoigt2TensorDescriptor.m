classdef SecondOrderVoigt2TensorDescriptor < handle
    
    properties (Access = protected, Abstract) 
        tensor
        indexTransformer
        dimVoigt
    end
    
    
    methods (Access = protected)
        
        function  representVoigtInTensor(obj,a)
            t = obj.tensor;
            converter = obj.indexTransformer;
            for iv = 1:obj.dimVoigt
                    [i,j] = converter.voigt2tensor(iv);
                    obj.voigtIndex = iv;
                    vf = obj.computeVoigtFactor();
                    aij = 1/vf*a(iv);                    
                    t(i,j) = aij;
                    t(j,i) = aij;

            end
            obj.tensor = t;
        end 
        
    end
        
    
    methods (Access = protected, Abstract)
        computeVoigtFactor(obj,index)
    end
    
    
    
    
    
end