classdef MonitoringFactory < handle

    methods (Access = public, Static)
        function m = create(cParams)
            switch cParams.type
                case 'OptimizationProblem'
                    m = MonitoringOptimizationProblem(cParams);
                case 'Null'
                    m = MonitoringNull(cParams);
            end
        end
    end
end