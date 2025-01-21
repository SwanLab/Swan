classdef testIsotropicFourthOrderTensor < handle
    
    properties (Access = public)
         tol = 1e-12;
    end
    
    properties (Access = private)
        E = 1;
        nu = 1/3;
        cVoigt
        cVoigtExplicit
    end
    
    methods (Access = public)
        
        function obj = testIsotropicFourthOrderTensor()
            obj.computeCheckedTensor();
            obj.computeToCheckTensor()
        end

        function error = computeError(obj)
            cV = obj.cVoigt.getValue();
            cVe = obj.cVoigtExplicit;
            error = norm(cV(:) - cVe(:));
        end

    end
    
    methods (Access = private)
        
        function computeCheckedTensor(obj)
          coef = obj.E/((obj.nu + 1)*(1 - 2*obj.nu));
          c1 = 1 - obj.nu;
          c2 = (1 - 2*obj.nu)/2;
          c3 = obj.nu;          
          C = [ c1  c3  c3   0   0   0;
                     c3  c1  c3   0   0   0;
                     c3  c3  c1   0   0   0;
                      0   0   0  c2   0   0;
                      0   0   0   0  c2   0;
                      0   0   0   0   0  c2];
          C = coef*C;
          obj.cVoigtExplicit = C;
        end
        
        function computeToCheckTensor(obj)
            C = IsotropicConstitutiveTensor(obj.E,obj.nu);
            obj.cVoigt = Tensor2VoigtConverter.convert(C);
        end

    end

end