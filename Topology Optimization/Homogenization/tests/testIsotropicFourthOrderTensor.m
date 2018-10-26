classdef testIsotropicFourthOrderTensor < test
    
    properties (Access = private)
        E = 1;
        nu = 1/3;
        ToCheckTensor
        CheckedTensor
    end
    
    methods (Access = public)
        
        function obj = testIsotropicFourthOrderTensor()
            obj.computeCheckedTensor();
            obj.computeToCheckTensor()
        end
    end
    
    methods (Access = private)
        
        function computeCheckedTensor(obj)        
          coef = obj.E/((obj.nu + 1)*(1 - 2*obj.nu));
          c1 = 1 - obj.nu;
          c2 = (1 - 2*obj.nu)/2;
          c3 = obj.nu;          
          Tensor = [ c1  c3  c3   0   0   0;
                     c3  c1  c3   0   0   0;
                     c3  c3  c1   0   0   0;
                      0   0   0  c2   0   0;
                      0   0   0   0  c2   0;
                      0   0   0   0   0  c2];          
          Tensor = coef*Tensor;
          obj.CheckedTensor = Tensor;          
        end
        
        function computeToCheckTensor(obj)
            Tensor = IsotropicConstitutiveTensor3D(obj.E,obj.nu);
            obj.ToCheckTensor = Tensor.tensorVoigt;
        end

    end
    
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            hasPassed = norm(double(obj.ToCheckTensor(:)) - obj.CheckedTensor(:)) < 1e-6;
        end
    end
end

