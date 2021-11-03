classdef testNotCommutingHomogPlaneStress < ...
            testCommutingHomogPlaneStress

    methods (Access = protected, Static)

        function nu = createPoissonValue()
            nu = 1/3;
        end
    end

    methods (Access = public)

        function hasPassed = hasPassed(obj)
            c1 = obj.vhpTensor.getValue();
            c2 = obj.vphTensor.getValue();
            error = norm(c2-c1)/norm(c1);
            hasPassed = error > 1e-6;
        end

    end

end