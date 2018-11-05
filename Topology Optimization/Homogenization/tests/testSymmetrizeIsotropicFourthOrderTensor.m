classdef testSymmetrizeIsotropicFourthOrderTensor < test
    
    properties (Access = protected)
          Ciso
          Csym
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
            obj.Csym = FourthOrderTensor();
            obj.Csym.setValue(obj.Ciso.getValue());
            obj.Csym.MakeMajorAndMinorSymmetrization();            
        end
        
      end
    
    methods (Access = protected)      
        function hasPassed = hasPassed(obj)
            Cs = obj.Csym.getValue();
            Ci = obj.Ciso.getValue();
            hasPassed = norm(Cs(:) - Ci(:)) < 1e-6;
        end
        
    end
end

