classdef testCommutingHomogPlaneStressWithZeroPoisson < ...
             testCommutingHomogPlaneStress

    properties (Access = public)
        tol = 1e-12;
    end

    methods (Access = protected, Static)
        
        function nu = createPoissonValue()
            nu = 0;
        end
    end
    
    methods (Access = public)
               
        function error = computeError(obj)
            c1 = obj.vhpTensor.getValue();
            c2 = obj.vphTensor.getValue();
            error = norm(c2-c1)/norm(c1);
        end
        
    end
    
end