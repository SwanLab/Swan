classdef MonitoringOptimizationProblem < handle

    properties (Access = private)
        cost
        constraint
        designVariable
        dualVariable
        problemParameters
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
            obj.initProblemParameters(cParams);
            obj.initOptimizationParameters(cParams);
        end

        function initProblemParameters(obj,cParams)

        end

        function initOptimizationParameters(obj,cParams)

        end

        function setNumberOfRowsAndColumns(obj)
            nF             = obj.cost.obtainNumberFields();
            nConstr        = obj.constraint.obtainNumberFields();
            nPlotsStandard = 2 + nF + 2*nConstr;
            nPlots         = nPlotsStandard; % +...
            obj.nRow       = floor(nPlots/7)+1;
            obj.nColumn    = min(nPlots,7);
        end

        function createMonitoring(obj)
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
            obj.createMonitoringOfVariable(2+i+j,'Norm L2 x',obj.designVariable)
        end

        function createMonitoringOfVariable(obj,i,title,f)
            chartType = obj.getChartType(title);
            newFig = DisplayFactory.create(chartType,title);
            obj.appendFigure(newFig);
            obj.figures{i}.show(obj.nRow,obj.nColumn,i,[0.06 0.06]);
            drawnow
            obj.data{i} = f;
        end

        function appendFigure(obj,fig)
            obj.figures{end+1} = fig;
        end

        function plot(obj,it)
            % cost
            % functionals along weights
            % diff constraints
            % l2 norm x
            % diff Lagrange multipliers
            % problemType plot
            % optType plot

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