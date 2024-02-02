classdef MonitoringOptimizationProblem < handle

    properties (Access = private)
        standardMonitoring
        functionalsMonitoring
        optParamsMonitoring
        maxNColumns
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
            obj.standardMonitoring    = MonitoringOptimizationProblemStandard(cParams);
            nStd                      = obj.standardMonitoring.nPlots;
            cParams.currentFigHandle  = nStd;
            obj.functionalsMonitoring = MonitoringOptimizationProblemFunctionals(cParams);
            nFun                      = obj.functionalsMonitoring.nPlots;
            cParams.currentFigHandle  = nStd+nFun;
            obj.optParamsMonitoring   = MonitoringOptimizationProblemParameters(cParams);
            obj.maxNColumns           = cParams.maxNColumns;
        end

        function setNumberOfRowsAndColumns(obj)
            nStd        = obj.standardMonitoring.nPlots;
            nFun        = obj.functionalsMonitoring.nPlots;
            nOpt        = obj.optParamsMonitoring.nPlots;
            nPlots      = nStd+nFun+nOpt;
            maxC        = obj.maxNColumns;
            obj.nRow    = ceil(nPlots/maxC);
            obj.nColumn = min(nPlots,maxC);
        end

        function createMonitoring(obj)
            m.figures = obj.figures;
            m.data    = obj.data;
            m.nRow    = obj.nRow;
            m.nColumn = obj.nColumn;

            figure
            m = obj.standardMonitoring.create(m);
            m = obj.functionalsMonitoring.create(m);
            m = obj.optParamsMonitoring.create(m);

            obj.figures = m.figures;
            obj.data    = m.data;
        end

        function plot(obj,it)
            m.figures = obj.figures;
            m.data    = obj.data;

            m = obj.standardMonitoring.plot(m,it);
            m = obj.functionalsMonitoring.plot(m,it);
            m = obj.optParamsMonitoring.plot(m,it);

            obj.figures = m.figures;
            obj.data    = m.data;
        end
    end
end