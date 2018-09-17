classdef SecondOrderTensor < handle
    
    properties
        tensor
        tensorVoigt
        tensorVoigtInPlaneStress
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
        
    end
    
    methods (Access = private)
        
        function generateRandomTensor(obj)
            obj.tensor = rand(3,3);
            obj.tensor = 0.5*(obj.tensor + obj.tensor');
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
    end
    
end

