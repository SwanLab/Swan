classdef testSymmetrizeIsotropicFourthOrderTensor < test
    
    properties (Access = protected)
          Ciso
          Csym
    end
    
    methods
        
        function obj = testSymmetrizeIsotropicFourthOrderTensor()
            obj.computeIsotropicFourthOrderTensor();
            obj.creatreSymmetricFourthOrderTensor()
        end
        
       
        
        function computeIsotropicFourthOrderTensor(obj)
            E = 1; nu = 1/3;
            obj.Ciso = IsotropicConstitutiveTensor3D(E,nu);
        end
        
        function creatreSymmetricFourthOrderTensor(obj)
            obj.Csym = fourthOrderTensor();
            obj.Csym.tensor = obj.Ciso.tensor;
            obj.Csym.MakeMajorAndMinorSymmetrization();            
        end
        
      end
    
    methods (Access = protected)      
        function hasPassed = hasPassed(obj)
            hasPassed = norm(obj.Csym.tensor(:) - obj.Ciso.tensor(:)) < 1e-6;
        end
        
    end
end

