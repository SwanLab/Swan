classdef FourthOrderVoigt2TensorConverter < Voigt2TensorConverterFor3DTensors
    

    methods (Access = public)
        
        function obj = FourthOrderVoigt2TensorConverter(tensor)
            obj.computeConversion(tensor)
        end        
    end
    
    methods (Access = protected)
        
        function selectTensorClass(obj)
            if obj.voigtTensor.getElasticityCase == '3D'
                obj.tensor = SymmetricFourthOrder3DTensor();
            elseif obj.voigtTensor.getElasticityCase == 'PlaneStress'
                obj.tensor = SymmetricFourthOrderPlaneStressTensor();
            end
        end
        
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
            obj.tensor.setValue(c);
        end
        
    end

    
    methods (Access = private,Static)
       
        function f = getVoigtFactor(iv,jv)
            if ((iv > 3 && jv <= 3) )  || ((iv <= 3 && jv > 3) )
                f = 1/2;
            elseif (iv>3 && jv>3)
                f = 1/4;
            else
                f = 1;
            end
        end

    end
end