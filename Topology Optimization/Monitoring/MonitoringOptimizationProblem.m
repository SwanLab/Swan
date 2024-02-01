classdef MonitoringOptimizationProblem < handle

    properties (Access = private)
        cost
        constraint
        designVariable
        dualVariable
        problemFunctionals
        optimizationParameters
    end

    properties (Access = private)
        standardMonitoring
        nRow
        nColumn
        figures
        data
    end

    methods (Access = public)
        function obj = MonitoringOptimizationProblem(cParams)
            obj.init(cParams);
            obj.standardMonitoring = MonitoringOptimizationProblemStandard(cParams);
            obj.setNumberOfRowsAndColumns();
            obj.createMonitoring();
        end

        function update(obj,it)
            obj.plot(it);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.designVariable = cParams.designVariable;
            obj.dualVariable   = cParams.dualVariable;
            obj.initProblemFunctionals(cParams);
            obj.initOptimizationParameters(cParams);
        end

        function initProblemFunctionals(obj,cParams)
            switch cParams.problemType
                case 'Topology'
                    obj.problemFunctionals{1} = cParams.volume;
            end
        end

        function initOptimizationParameters(obj,cParams)
            switch cParams.optimizationType
                case 'NullSpace'
                    obj.optimizationParameters.type   = 'NullSpace';
                    obj.optimizationParameters.primal = cParams.primalUpdater;
                    obj.optimizationParameters.nPlots = 2;
                otherwise
                    obj.optimizationParameters.type   = 'Null';
                    obj.optimizationParameters.nPlots = 0;
            end
        end

        function setNumberOfRowsAndColumns(obj)
            nPlotsStd      = obj.standardMonitoring.nPlots;
            nPlotsProblem  = length(obj.problemFunctionals);
            nPlotsOpt      = obj.optimizationParameters.nPlots;
            nPlots         = nPlotsStd+nPlotsProblem+nPlotsOpt;
            obj.nRow       = ceil(nPlots/7);
            obj.nColumn    = min(nPlots,7);
        end

        function createMonitoring(obj)
            m.figures = obj.figures;
            m.data    = obj.data;
            m.nRow    = obj.nRow;
            m.nColumn = obj.nColumn;
            figure
            m = obj.standardMonitoring.create(m);

            obj.figures = m.figures;
            obj.data    = m.data;

            obj.createProblemFunctionalsMonitoring(); % another class
            obj.createOptimizationParametersMonitoring(); % another class
        end

        function createProblemFunctionalsMonitoring(obj)
            n = length(obj.data);
            for i = 1:length(obj.problemFunctionals)
                J     = obj.problemFunctionals{i};
                title = J.getTitleToPlot();
                obj.createMonitoringOfVariable(n+i,title,J);
            end
        end

        function createOptimizationParametersMonitoring(obj)
            n = length(obj.data);
            switch obj.optimizationParameters.type
                case 'Null'
                case 'NullSpace'
                    primal = obj.optimizationParameters.primal;
                    obj.createMonitoringOfVariable(n+1,'Line Search',primal);
                    obj.createMonitoringOfVariable(n+2,'Line Search trials',primal);
            end
        end

        function createMonitoringOfVariable(obj,i,title,f)
            chartType = obj.getChartType(title);
            newFig = DisplayFactory.create(chartType,title);
            obj.appendFigure(newFig);
            obj.figures{i}.show(obj.nRow,obj.nColumn,i,[0.06 0.04]);
            drawnow
            obj.data{i} = f;
        end

        function appendFigure(obj,fig)
            obj.figures{end+1} = fig;
        end

        function plot(obj,it)
            m.figures = obj.figures;
            m.data    = obj.data;
            obj.standardMonitoring.plot(m,it);

            obj.figures = m.figures;
            obj.data    = m.data;

            obj.plotProblemFunctionals(it); %
            obj.plotOptimizationParameters(it); %
        end

        function plotProblemFunctionals(obj,it)
            x       = obj.designVariable;
            nF      = obj.cost.obtainNumberFields();
            nConstr = obj.constraint.obtainNumberFields();
            n       = 2+nF+2*nConstr;
            for i = 1:length(obj.problemFunctionals)
                [J,~] = obj.data{n+i}.computeFunctionAndGradient(x);
                obj.figures{n+i}.updateParams(it,J);
                obj.figures{n+i}.refresh();
            end
        end

        function plotOptimizationParameters(obj,it)
            nF      = obj.cost.obtainNumberFields();
            nConstr = obj.constraint.obtainNumberFields();
            nFunct  = length(obj.problemFunctionals);
            n       = 2+nF+2*nConstr+nFunct;
            switch obj.optimizationParameters.type
                case 'Null'
                case 'NullSpace'
                    obj.figures{n+1}.updateParams(it,obj.data{n+1}.tau);
                    obj.figures{n+2}.updateParams(it,obj.data{n+2}.getCurrentTrials());
                    obj.figures{n+1}.refresh();
                    obj.figures{n+2}.refresh();
            end
        end
    end

    methods (Static, Access = private)
        function type = getChartType(title)
            switch title
                case {'Line Search','Line Search trials'}
                    type = 'bar';
                case {'Norm L2 x'}
                    type = 'log';
                otherwise 
                    type = 'plot';
            end
        end
    end
end