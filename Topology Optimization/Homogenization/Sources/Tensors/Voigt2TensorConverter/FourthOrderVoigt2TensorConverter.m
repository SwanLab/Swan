classdef FourthOrderVoigt2TensorConverter < Voigt2TensorConverterFor3DTensors
    

    methods (Access = public)
        
        function obj = FourthOrderVoigt2TensorConverter(tensor)
            obj.computeConversion(tensor)
        end        
    end
    
    methods (Access = protected)
        
        function  representVoigtInTensor(obj,a)
            c = obj.tensor;
            converter = obj.indexTransformer;
            for iv = 1:obj.dimVoigt
                for jv = 1:obj.dimVoigt
                    [i,j] = converter.voigt2tensor(iv);
                    [k,l] = converter.voigt2tensor(jv);
                    
                    if ((iv > 3 && jv <= 3) )  || ((iv <= 3 && jv > 3) )
                        %factor =1/sqrt(2);
                        factor =1/2;
                    elseif (iv>3 && jv>3) 
                        factor = 1/4;
                    else
                        factor = 1;
                    end
                    aij = factor*a(iv,jv);
                    
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
            obj.tensor = c;
        end
        
        
        function obtainTensorSize(obj)
            obj.tensorSize = [obj.dim,obj.dim,obj.dim,obj.dim];
        end
        
    end
    
end