classdef SecondOrderVoigt2TensorConverter < Voigt2TensorConverterFor3DTensors
    
    properties (Access = protected)
        voigtIndex
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
                
        function obtainTensorSize(obj)
            obj.tensorSize = [obj.dim,obj.dim];
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeVoigtFactor(obj,index)
    end
    
end