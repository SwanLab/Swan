classdef AnisotropicContributionTensor < FourthOrderTensor
    
    properties (Access = public)
        FirstTensor
        SecondTensor
        ThirdTensor
    end
    
    methods (Access = public)
        
        function obj = AnisotropicContributionTensor(A,direction)
            obj.generateTensors(A,direction);
            obj.addTensors()
        end
    end
    
    methods (Access = private)
        
        function generateTensors(obj,A,direction)
            obj.FirstTensor  = A;
            obj.SecondTensor = SecondComplementaryTensor(A,direction);
            obj.ThirdTensor  = ThirdComplementaryTensor(A,direction);
        end

        function addTensors(obj)
            t1  = obj.FirstTensor.tensor;
            t2  = obj.SecondTensor.tensor;
            t3  = obj.ThirdTensor.tensor;
            obj.tensor = t1 + t2 + t3;
        end
        
    end
    
end

