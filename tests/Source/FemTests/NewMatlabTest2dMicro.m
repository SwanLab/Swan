classdef NewMatlabTest2dMicro < NewMatlabTest

    methods (Access = protected)
        function computeVariableThroughFemSolver(obj)
            femSolver = Elastic_Problem_Micro.create(obj.testName);
            props.kappa = .9107;
            props.mu    = .3446;
            femSolver.setMatProps(props);
            femSolver.computeChomog;
            obj.fem = femSolver;
        end        
    end

end

