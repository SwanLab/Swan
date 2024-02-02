classdef MonitoringOptimizationProblemFunctionals < handle

    properties (Access = public)
        nPlots
    end

    properties (Access = private)
        designVariable
        problemFunctionals
        currentFigHandle
    end

    methods (Access = public)
        function obj = MonitoringOptimizationProblemFunctionals(cParams)
            obj.init(cParams);
        end

        function m = create(obj,m)
            n = obj.currentFigHandle;
            for i = 1:length(obj.problemFunctionals)
                J     = obj.problemFunctionals{i};
                title = J.getTitleToPlot();
                m = MonitoringVariable.create(m,n+i,title,J);
            end
        end

        function m = plot(obj,m,it)
            x = obj.designVariable;
            n = obj.currentFigHandle;
            for i = 1:length(obj.problemFunctionals)
                [J,~] = m.data{n+i}.computeFunctionAndGradient(x);
                m.figures{n+i}.updateParams(it,J);
                m.figures{n+i}.refresh();
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.designVariable     = cParams.designVariable;
            obj.problemFunctionals = cParams.functionals;
            obj.nPlots             = length(obj.problemFunctionals);
            obj.currentFigHandle   = cParams.currentFigHandle;
        end
    end
end