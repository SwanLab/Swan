classdef FourthOrderTensor2VoigtConverter < Tensor2VoigtConverterFor3DTensors
    

    methods (Access = public)
        
        function obj = FourthOrderTensor2VoigtConverter(Tensor)
            obj.computeConversion(Tensor)
        end        
    end
    
    methods (Access = protected)
                
       function  representTensorInVoigt(obj,A)
            c = obj.voigtTensor;
            converter = obj.indexTransformer; 
            for i = 1:obj.dim
                for j = 1:obj.dim
                    for k = 1:obj.dim
                        for l = 1:obj.dim
                            iv = converter.tensor2Voigt(i,j);
                            jv = converter.tensor2Voigt(k,l);
                            
                            
                            if ((iv > 3 && jv <= 3) )  || ((iv <= 3 && jv > 3) )
                                %factor =sqrt(2);
                                factor =1;
                            elseif (iv>3 && jv>3)
                                %factor = 2;
                                factor = 1;
                            else
                                factor = 1;
                            end
                            
                            
                            c(iv,jv) = factor*A(i,j,k,l);
                        end
                    end
                end
            end
            obj.voigtTensor = c;
       end
        
       
       function obtainVoigtTensorSize(obj)
            obj.voigtTensorSize = [obj.dimVoigt,obj.dimVoigt]; 
       end
       
    end
    
end

