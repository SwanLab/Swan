classdef SecondOrderTensor < handle
    
    properties (Access = public)
        tensor
        tensorVoigt
        tensorVoigtInPlaneStress
    end
    
    properties (Abstract)
        VoigtFactor
    end
    
    methods
        
        function obj = SecondOrderTensor()
            obj.generateRandomTensor();
            obj.transformTensor2Voigt();
        end
        
        
        function  transformTensor2Voigt(obj)
            obj.tensorVoigt = zeros(6,1);
            for i = 1:size(obj.tensor,1)
                for j = 1:size(obj.tensor,2)
                    iv = obj.transformTensorIndex2VoigtIndex(i,j);
                    obj.computeVoigtFactor(iv);
                    obj.tensorVoigt(iv) = obj.VoigtFactor*obj.tensor(i,j);
                end
            end
        end
        
        function transformVoigt2Tensor(obj)
            obj.tensor = zeros(3,3);
            for iVoigt = 1:length(obj.tensorVoigt)
                    [iv,jv] = obj.transformVoigtIndex2TensorIndex(iVoigt);
                    obj.computeVoigtFactor(iVoigt);
                    obj.tensor(iv,jv) = 1/obj.VoigtFactor*obj.tensorVoigt(iVoigt);
                    obj.tensor(jv,iv) = 1/obj.VoigtFactor*obj.tensorVoigt(iVoigt);                    
            end
        end
    end
    

    methods (Access = private)
        
        function generateRandomTensor(obj)
            obj.tensor = rand(3,3);
            obj.tensor = 0.5*(obj.tensor + obj.tensor');
        end
        
    end
    
    
    methods (Access = protected)
        
        function computePlaneStressTensor(obj)
            obj.tensorVoigtInPlaneStress(1,1) = obj.tensorVoigt(1);
            obj.tensorVoigtInPlaneStress(2,1) = obj.tensorVoigt(2);
            obj.tensorVoigtInPlaneStress(3,1) = obj.tensorVoigt(6);
        end
    end
    
    methods (Abstract)
        computeVoigtFactor(obj,a)
    end
    

    
    methods (Static)
        function ind = transformTensorIndex2VoigtIndex(istre,jstre)
            Tensor2Voigt =  [1 6 5;
                             6 2 4;
                             5 4 3];
            
            ind = Tensor2Voigt(istre,jstre);
        end
   
        
        function [istre,jstre] = transformVoigtIndex2TensorIndex(VoigtIndex)
            Voigt2Tensor =  [1 1;
                             2 2;
                             3 3
                             2 3;
                             1 3;
                             1 2];
                         
            istre = Voigt2Tensor(VoigtIndex,1);             
            jstre = Voigt2Tensor(VoigtIndex,2);             
            
        end
 
        
        
        
        
        function isSymmetric = isSymmetric(A)
            Asym = SecondOrderTensor.symmetrize(A);
            isSymmetric = norm(Asym(:)) - norm(A(:)) < 1e-14;
        end


        
        function Asym = symmetrize(A)
            Asym = 0.5*(A+A');
        end
        
        
        
    end
    
end

