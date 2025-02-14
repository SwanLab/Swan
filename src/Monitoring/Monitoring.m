classdef Monitoring < handle

    properties (Access = private)
        figures
    end

    properties (Access = private)
        shallDisplay
        maxNColumns
        titles
        chartTypes
    end

    methods (Access = public)
        function obj = Monitoring(cParams)
            obj.init(cParams);
            obj.createMonitoring(cParams);
        end

        function update(obj,it,data)
            nPlots = length(obj.figures);
            for i = 1:nPlots
                obj.figures{i}.updateParams(it,data{i});
            end
        end

        function refresh(obj)
            nPlots = length(obj.figures);
            for i = 1:nPlots
                obj.figures{i}.refresh();
            end      
            drawnow;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.shallDisplay  = cParams.shallDisplay;
            obj.maxNColumns   = cParams.maxNColumns;
            obj.titles        = cParams.titles;
            obj.chartTypes    = cParams.chartTypes;
        end

        function [nRow,nColumn] = computeNumberRowsColumns(obj)
            nPlots  = length(obj.titles);
            maxC    = obj.maxNColumns;
            nRow    = ceil(nPlots/maxC);
            nColumn = min(nPlots,maxC);
        end

        function createMonitoring(obj,cParams)
            if (obj.shallDisplay)
                figure
                nPlots         = length(obj.titles);
                [nRow,nColumn] = obj.computeNumberRowsColumns();
                idxMultiBar = 1; idxSurf = 1;
                for i = 1:nPlots
                    sDisp.title     = obj.titles{i};
                    sDisp.chartType = obj.chartTypes{i};
                    sDisp.figureIdx = length(findobj('Type','Figure'));
                    if sDisp.chartType == "multiplot"
                        sDisp.legend = cParams.legends{idxMultiBar};
                        idxMultiBar = idxMultiBar+1;
                    elseif sDisp.chartType == "surf"
                        sDisp.barLim = cParams.barLims{idxSurf};
                        sDisp.fun    = cParams.funs{idxSurf};
                        idxSurf = idxSurf+1;
                    end
                    sDisp.position  = i;
                    newFig    = DisplayAbstract.create(sDisp);
                    obj.appendFigure(newFig);
                    obj.figures{i}.show(nRow,nColumn,i,[0.06 0.04]);
                    hold on
                end
            end
        end

        function appendFigure(obj,fig)
            obj.figures{end+1} = fig;
        end

    end
end