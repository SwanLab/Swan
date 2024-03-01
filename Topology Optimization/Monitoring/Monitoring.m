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
            obj.createMonitoring();
        end

        function update(obj,it,data)
            if (obj.shallDisplay)
                obj.plot(it,data);
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.shallDisplay = cParams.shallDisplay;
            obj.maxNColumns  = cParams.maxNColumns;
            obj.titles       = cParams.titles;
            obj.chartTypes   = cParams.chartTypes;
        end

        function [nRow,nColumn] = computeNumberRowsColumns(obj)
            nPlots  = length(obj.titles);
            maxC    = obj.maxNColumns;
            nRow    = ceil(nPlots/maxC);
            nColumn = min(nPlots,maxC);
        end

        function createMonitoring(obj)
            if (obj.shallDisplay)
                figure
                nPlots         = length(obj.titles);
                [nRow,nColumn] = obj.computeNumberRowsColumns();
                for i = 1:nPlots
                    title     = obj.titles{i};
                    chartType = obj.chartTypes{i};
                    newFig    = DisplayFactory.create(chartType,title);
                    obj.appendFigure(newFig);
                    obj.figures{i}.show(nRow,nColumn,i,[0.06 0.04]);
                    hold on
                end
            end
        end

        function appendFigure(obj,fig)
            obj.figures{end+1} = fig;
        end

        function plot(obj,it,data)
            nPlots = length(obj.figures);
            for i = 1:nPlots
                obj.figures{i}.updateParams(it,data(i));
                obj.figures{i}.refresh();
            end
        end
    end
end