classdef testFourthOrderTensor < test
    
    properties (Access = protected)
        E = 1;
        nu = 1/3;
        ToCheckTensor
        CheckedTensor
    end
    
    methods
        
        function obj = testFourthOrderTensor()
            obj.computeCheckedTensor();
            obj.computeToCheckTensor()
        end
        
        function computeCheckedTensor(obj)
          Tensor = zeros(6,6);          
          coef = obj.E/((obj.nu + 1)*(1 - 2*obj.nu));
          coef1 = 1 - obj.nu;
          coef2 = (1 - 2*obj.nu)/2;
          
          Tensor(1,1) = coef1;
          Tensor(2,2) = coef1;
          Tensor(3,3) = coef1;
          
          Tensor(4,4) = coef2;
          Tensor(5,5) = coef2;
          Tensor(6,6) = coef2;
          
          Tensor(1,2) = obj.nu;
          Tensor(2,1) = obj.nu;
          Tensor(1,3) = obj.nu;
          Tensor(3,1) = obj.nu;
          Tensor(3,2) = obj.nu;
          Tensor(2,3) = obj.nu;
          
          Tensor = coef*Tensor;
          obj.CheckedTensor = Tensor;
          
        end
        
        function computeToCheckTensor(obj)
            Tensor = IsotropicConstitutiveTensor3D(obj.E,obj.nu);
            %Tensor.tensorVoigt = zeros(6,6);
            obj.ToCheckTensor = Tensor.tensorVoigt;
        end

    end
    
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            hasPassed = norm(double(obj.ToCheckTensor(:)) - obj.CheckedTensor(:)) < 1e-6;
        end
    end
end

