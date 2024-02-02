classdef MonitoringOptimizationProblemParameters < handle

    properties (Access = public)
        nPlots
    end

    properties (Access = private)
        type
        parameters
        currentFigHandle
    end

    methods (Access = public)
        function obj = MonitoringOptimizationProblemParameters(cParams)
            obj.init(cParams);
        end

        function m = create(obj,m)
            n = obj.currentFigHandle;
            switch obj.type
                case 'NullSpace'
                    primal = obj.parameters.primal;
                    m = MonitoringVariable.create(m,n+1,'Line Search',primal);
                    m = MonitoringVariable.create(m,n+2,'Line Search trials',primal);
                %case...
            end
        end

        function m = plot(obj,m,it)
            n = obj.currentFigHandle;
            switch obj.type
                case 'NullSpace'
                    m.figures{n+1}.updateParams(it,m.data{n+1}.tau);
                    m.figures{n+2}.updateParams(it,m.data{n+2}.getCurrentTrials());
                    m.figures{n+1}.refresh();
                    m.figures{n+2}.refresh();
                %case...
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.currentFigHandle = cParams.currentFigHandle;
            obj.type             = cParams.optimizationType;
            switch obj.type
                case 'NullSpace'
                    obj.parameters.primal = cParams.primalUpdater;
                    obj.nPlots = 2;
                %case...
                otherwise
                    obj.nPlots = 0;
            end
        end
    end
end