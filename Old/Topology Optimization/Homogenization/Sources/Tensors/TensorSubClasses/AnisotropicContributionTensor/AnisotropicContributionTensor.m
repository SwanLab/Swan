classdef AnisotropicContributionTensor < Stiffness3DTensor
    
    properties (Access = private)
        firstTensor
        secondTensor
        thirdTensor
        dir
    end
    
    methods (Access = public)
        
        function obj = AnisotropicContributionTensor(A,direction)
            obj.generateTensors(A,direction);
            obj.addTensors()
        end
        
        function t = clone(obj)
            tens = obj.firstTensor();
            d = obj.dir;
            t = AnisotropicContributionTensor(tens,d);
        end
    end
    
    methods (Access = private)
        
        function generateTensors(obj,A,direction)
            obj.dir          = direction;
            obj.firstTensor  = A;
            obj.secondTensor = SecondAnisotropicTensor(A,obj.dir);
            obj.thirdTensor  = ThirdAnisotropicTensor(A,obj.dir);
        end

        function addTensors(obj)
            t1  = obj.firstTensor.getValue();
            t2  = obj.secondTensor.getValue();
            t3  = obj.thirdTensor.getValue();
            tensor = t1 + t2 + t3;
            obj.setValue(tensor);
        end
        
    end
    
end