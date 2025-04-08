classdef MonitoringProjectedGradient < handle

    properties (Access = private)
        cost
        designVariable
        primalUpdater
    end

    properties (Access = private)
        monitoring
    end

    methods (Access = public)
        function obj = MonitoringProjectedGradient(cParams)
            obj.init(cParams);
            obj.createMonitoring(cParams);
        end

        function update(obj,nIter,sD)
            data = obj.computeDataUpdated(nIter,sD);
            obj.monitoring.update(nIter,num2cell(data));
        end

        function refresh(obj)
            obj.monitoring.refresh();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.cost             = cParams.cost;
            obj.designVariable   = cParams.designVariable;
            obj.primalUpdater    = cParams.primalUpdater;
        end

        function createMonitoring(obj,cParams)
            titlesF       = obj.cost.getTitleFields();
            nSFCost       = length(titlesF);
            titles        = [{'Cost'};titlesF;{'Norm L2 x'}];
            titles  = [titles;{'Line Search';'Line Search trials';'Merit'}];
            chCost = cell(1,nSFCost);
            for i = 1:nSFCost
                chCost{i} = 'plot';
            end
            chartTypes = [{'plot'},chCost,{'logy'},{'bar','bar','plot'}];
            s.shallDisplay = cParams.shallDisplay;
            s.maxNColumns  = 6;
            s.titles       = titles;
            s.chartTypes   = chartTypes;
            obj.monitoring = Monitoring(s);
        end

        function data = computeDataUpdated(obj,nIter,sD)
            data = obj.cost.value;
            data = [data;obj.cost.getFields(':')];
            data = [data;obj.designVariable.computeL2normIncrement()];
            if nIter == 0
                data = [data;0;0;0];
            else
                data = [data;obj.primalUpdater.tau;sD.lineSearchTrials;sD.meritNew];
            end
        end
    end
end