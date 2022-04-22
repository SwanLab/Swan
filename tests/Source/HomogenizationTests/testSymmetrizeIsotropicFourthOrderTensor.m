classdef testSymmetrizeIsotropicFourthOrderTensor < handle
    
    properties (Access = public)
          tol = 1e-6;
    end

    properties (Access = protected)
          Ciso
          Csym
    end
    
    methods (Access = public)
        
        function obj = testSymmetrizeIsotropicFourthOrderTensor()
            obj.computeIsotropicFourthOrderTensor();
            obj.createSymmetricFourthOrderTensor()
        end

        function error = computeError(obj)
            Cs = obj.Csym.getValue();
            Ci = obj.Ciso.getValue();
            error = norm(Cs(:) - Ci(:));
        end
    end

    methods (Access = private)

        function computeIsotropicFourthOrderTensor(obj)
            E = 1; nu = 1/3;
            obj.Ciso = IsotropicConstitutiveTensor(E,nu);
        end
        
        function createSymmetricFourthOrderTensor(obj)
            obj.Csym = SymmetricFourthOrder3DTensor();
            obj.Csym.setValue(obj.Ciso.getValue());
            obj.Csym.MakeMajorAndMinorSymmetrization();
        end

    end

end
