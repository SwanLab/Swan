classdef MonitoringOptimizationProblemStandard < handle

    properties (Access = public)
        nPlots
    end

    properties (Access = private)
        cost
        constraint
        designVariable
        dualVariable
        isConstrained
    end

    methods (Access = public)
        function obj = MonitoringOptimizationProblemStandard(cParams)
            obj.init(cParams);
        end

        function m = create(obj,m)
            m = obj.createMonitoring(m);
        end

        function m = plot(obj,m,it)
            nF      = obj.cost.obtainNumberFields();
            nConstr = obj.constraint.obtainNumberFields();
            m.figures{1}.updateParams(it,m.data{1}.value);
            m.figures{1}.refresh();
            for i = 1:nF
                m.figures{i+1}.updateParams(it,m.data{i+1}.getFields(i));
                m.figures{i+1}.refresh();
            end
            for j = 1:nConstr
                m.figures{j+i+1}.updateParams(it,m.data{j+i+1}.value(j,1));
                m.figures{j+i+1}.refresh();
            end
            m.figures{j+i+2}.updateParams(it,m.data{j+i+2}.computeL2normIncrement());
            m.figures{j+i+2}.refresh();
            for k = 1:nConstr
                m.figures{k+j+i+2}.updateParams(it,m.data{k+j+i+2}.value(k,1));
                m.figures{k+j+i+2}.refresh();
            end
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.designVariable = cParams.designVariable;
            obj.isConstrained  = cParams.isConstrained;
            switch cParams.isConstrained
                case true
                    nConstr          = cParams.constraint.obtainNumberFields();
                    obj.constraint   = cParams.constraint;
                    obj.dualVariable = cParams.dualVariable;
                case false
                    nConstr  = 0;
            end
            nF         = obj.cost.obtainNumberFields();
            obj.nPlots = 2+nF+2*nConstr;
        end

        function m = createMonitoring(obj,m)
            m = obj.monitorCost(m);
            m = obj.monitorConstraints(m);
            m = obj.monitorDesignVariable(m);
            m = obj.monitorLagrangeMultipliers(m);
        end

        function m = monitorCost(obj,m)
            titlesF = obj.cost.getTitleFields();
            m       = obj.createMonitoringOfVariable(m,1,'Cost',obj.cost);
            for i = 1:length(titlesF)
                m = obj.createMonitoringOfVariable(m,1+i,titlesF{i},obj.cost);
            end
        end

        function m = monitorConstraints(obj,m)
            if obj.isConstrained
                titlesConst = obj.constraint.getTitleFields();
                n           = length(m.figures);
                for i = 1:length(titlesConst)
                    m = obj.createMonitoringOfVariable(m,n+i,titlesConst{i},obj.constraint);
                end
            end
        end

        function m = monitorDesignVariable(obj,m)
            n = length(m.figures);
            m = obj.createMonitoringOfVariable(m,n+1,'Norm L2 x',obj.designVariable);
        end

        function m = monitorLagrangeMultipliers(obj,m)
            if obj.isConstrained
                titlesConst = obj.constraint.getTitleFields();
                n           = length(m.figures);
                for i = 1:length(titlesConst)
                    m = obj.createMonitoringOfVariable(m,n+1,['\lambda_{',titlesConst{i},'}'],obj.dualVariable);
                end
            end
        end

        function m = createMonitoringOfVariable(obj,m,i,title,f)
            chartType = obj.getChartType(title);
            newFig = DisplayFactory.create(chartType,title);
            m = obj.appendFigure(m,newFig);
            m.figures{i}.show(m.nRow,m.nColumn,i,[0.06 0.04]);
            drawnow
            m.data{i} = f;
        end
    end

    methods (Static, Access = private)
        function m = appendFigure(m,fig)
            m.figures{end+1} = fig;
        end

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