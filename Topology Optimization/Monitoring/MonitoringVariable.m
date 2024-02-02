classdef MonitoringVariable < handle

    methods (Static, Access = public)
        function m = create(m,iFig,title,varObj)
            obj = MonitoringVariable();
            chartType = obj.getChartType(title);
            newFig = DisplayFactory.create(chartType,title);
            m = obj.appendFigure(m,newFig);
            m.figures{iFig}.show(m.nRow,m.nColumn,iFig,[0.06 0.04]);
            drawnow
            m.data{iFig} = varObj;
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