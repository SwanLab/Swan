classdef testCommutingHomogPlaneStressWithZeroPoisson < ...
             testCommutingHomogPlaneStress &...
             testShowingError
    
    properties (Access = protected)
        tol = 1e-12;
    end
    
    methods (Access = protected, Static)
        
        function nu = createPoissonValue()
            nu = 0;
        end
    end
    
    methods (Access = protected)
               
        function computeError(obj)
            c1 = obj.vhpTensor.getValue();
            c2 = obj.vphTensor.getValue();
            obj.error = norm(c2-c1)/norm(c1);
        end
        
    end
    
end

