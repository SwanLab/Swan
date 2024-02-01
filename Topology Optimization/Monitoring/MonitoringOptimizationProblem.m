classdef MonitoringOptimizationProblem < handle  % Pending to divide in subClasses (std; functionals; optParams)

    properties (Access = private)
        cost
        constraint
        designVariable
        dualVariable
        problemFunctionals
        optimizationParameters
    end

    properties (Access = private)
        nRow
        nColumn
        figures
        data
    end

    methods (Access = public)
        function obj = MonitoringOptimizationProblem(cParams)
            obj.init(cParams);
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
            nF             = obj.cost.obtainNumberFields();
            nConstr        = obj.constraint.obtainNumberFields();
            nPlotsStandard = 2 + nF + 2*nConstr;
            nPlotsProblem  = length(obj.problemFunctionals);
            nPlotsOpt      = obj.optimizationParameters.nPlots;
            nPlots         = nPlotsStandard+nPlotsProblem+nPlotsOpt;
            obj.nRow       = ceil(nPlots/7);
            obj.nColumn    = min(nPlots,7);
        end

        function createMonitoring(obj)
            obj.createStandardMonitoring();
            obj.createProblemFunctionalsMonitoring();
            obj.createOptimizationParametersMonitoring();
        end

        function createStandardMonitoring(obj)
            titlesF     = obj.cost.getTitleFields();
            titlesConst = obj.constraint.getTitleFields();
            figure
            obj.createMonitoringOfVariable(1,'Cost',obj.cost);
            for i = 1:length(titlesF)
                obj.createMonitoringOfVariable(1+i,titlesF{i},obj.cost);
            end
            for j = 1:length(titlesConst)
                obj.createMonitoringOfVariable(1+i+j,titlesConst{j},obj.constraint);
            end
            obj.createMonitoringOfVariable(2+i+j,'Norm L2 x',obj.designVariable);
            for k = 1:length(titlesConst)
                obj.createMonitoringOfVariable(2+i+j+k,['\lambda_{',titlesConst{j},'}'],obj.dualVariable);
            end
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
            obj.plotStandard(it);
            obj.plotProblemFunctionals(it);
            obj.plotOptimizationParameters(it);
        end

        function plotStandard(obj,it)
            nF      = obj.cost.obtainNumberFields();
            nConstr = obj.constraint.obtainNumberFields();
            obj.figures{1}.updateParams(it,obj.data{1}.value);
            obj.figures{1}.refresh();
            for i = 1:nF
                obj.figures{i+1}.updateParams(it,obj.data{i+1}.getFields(i));
                obj.figures{i+1}.refresh();
            end
            for j = 1:nConstr
                obj.figures{j+i+1}.updateParams(it,obj.data{j+i+1}.value(j,1));
                obj.figures{j+i+1}.refresh();
            end
            obj.figures{j+i+2}.updateParams(it,obj.data{j+i+2}.computeL2normIncrement());
            obj.figures{j+i+2}.refresh();
            for k = 1:nConstr
                obj.figures{k+j+i+2}.updateParams(it,obj.data{k+j+i+2}.value(k,1));
                obj.figures{k+j+i+2}.refresh();
            end
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