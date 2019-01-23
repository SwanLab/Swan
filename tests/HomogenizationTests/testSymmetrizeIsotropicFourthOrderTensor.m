classdef testSymmetrizeIsotropicFourthOrderTensor < testShowingError
    
    properties (Access = protected)
          Ciso
          Csym
          tol = 1e-6;
    end
    
    methods
        
        function obj = testSymmetrizeIsotropicFourthOrderTensor()
            obj.computeIsotropicFourthOrderTensor();
            obj.createSymmetricFourthOrderTensor()
        end
        
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
    
    methods (Access = protected)      
        function computeError(obj)
            Cs = obj.Csym.getValue();
            Ci = obj.Ciso.getValue();
            obj.error = norm(Cs(:) - Ci(:));;
        end
        
    end
end
