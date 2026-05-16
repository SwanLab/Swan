classdef MonitoringNullSpace < handle

    properties (Access = private)
        cost
        constraint
        designVariable
        dualVariable
        primalUpdater
    end

    properties (Access = private)
        monitoring
    end

    methods (Access = public)
        function obj = MonitoringNullSpace(cParams)
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
            obj.constraint       = cParams.constraint;
            obj.designVariable   = cParams.designVariable;
            obj.dualVariable     = cParams.dualVariable;
            obj.primalUpdater    = cParams.primalUpdater;
        end

        function createMonitoring(obj,cParams)
            titlesF       = obj.cost.getTitleFields();
            %titlesF = [titlesF; {'Cost Input'}; {'Cost Output'}]; % add the input/output cost title
            
            titlesConst   = obj.constraint.getTitleFields();
            nSFCost       = length(titlesF);
            nSFConstraint = length(titlesConst);
            titles        = [{'Cost'};titlesF;titlesConst;{'Norm L2 x'}];
            chConstr      = cell(1,nSFConstraint);
            for i = 1:nSFConstraint
                titles{end+1} = ['\lambda_{',titlesConst{i},'}'];
                chConstr{i}   = 'plot';
            end
            titles  = [titles;{'Line Search';'Line Search trials';'Eta';'EtaMax';'lG';'lJ';'Merit'}];
            chCost = cell(1,nSFCost);
            for i = 1:nSFCost
                chCost{i} = 'plot';
            end
            chartTypes = [{'plot'},chCost,chConstr,{'logy'},chConstr,{'bar','bar','plot','logy','plot','plot','plot'}];
            switch class(obj.designVariable)
                case 'LevelSet'
                    titles = [titles;{'Theta';'Alpha';'Beta'}];
                    chartTypes = [chartTypes,{'plot','plot','plot'}];
            end
            s.shallDisplay = cParams.shallDisplay;
            s.maxNColumns  = 6;
            s.titles       = titles;
            s.chartTypes   = chartTypes;
            obj.monitoring = Monitoring(s);
        end

        function data = computeDataUpdated(obj,nIter,sD)
            data = obj.cost.value;
            data = [data;obj.cost.getFields(':')];

            % c_in  = 0;
            % c_out = 0;
            % 
            % shFun = obj.cost.shapeFunctions;
            % if iscell(shFun)
            %     for i = 1:length(shFun)
            %         if isa(shFun{i}, 'NonSelfAdjointComplianceFunctional')
            %             c_in  = shFun{i}.cost_in_norm;
            %             c_out = shFun{i}.cost_out_norm;
            %             break;
            %         end
            %     end
            % end
            % 
            % data = [data; c_in; c_out]; % add the data we want to plot
            data = [data;obj.constraint.value];
            data = [data;obj.designVariable.computeL2normIncrement()];
            data = [data;obj.dualVariable.fun.fValues];
            if nIter == 0
                data = [data;0;0;0;sD.etaMax;0;0;NaN];
            else
                data = [data;obj.primalUpdater.tau;sD.lineSearchTrials;sD.eta;sD.etaMax;norm(sD.lG);norm(sD.lJ);sD.meritNew];
            end
            switch class(obj.designVariable)
                case 'LevelSet'
                    if nIter == 0
                        data = [data;0;0;0];
                    else
                        data = [data;obj.primalUpdater.Theta;obj.primalUpdater.Alpha;obj.primalUpdater.Beta];
                    end
            end
        end
    end
end