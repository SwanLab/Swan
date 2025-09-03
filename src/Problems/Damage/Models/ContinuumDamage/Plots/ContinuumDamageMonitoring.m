classdef ContinuumDamageMonitoring < handle

    properties (GetAccess = public, SetAccess = private)
        data
    end

    properties (Access = private)
        monitor
        printInfo
    end

    methods (Access = public)

        function obj = ContinuumDamageMonitoring(cParams)
            obj.init(cParams);
        end
    end

    methods (Access = public)

        function update(obj,iter,data)
            obj.monitor.update(iter,data)
        end

        function refresh(obj)
            obj.monitor.refresh;
        end

        function updateAndRefresh(obj,iter,data)
            obj.update(iter,data);
            obj.refresh();
        end

        function printCost(obj,name,iter,cost,err)
            if obj.printInfo == true
                X = sprintf('%s:%d / cost: %.8e  (diff:%.8e) \n',name,iter,cost,err);
                fprintf(X);
            end
        end

        function printStep(obj,step,maxSteps)
            if obj.printInfo == true
                fprintf('\n ********* STEP %i/%i *********  \n',step,maxSteps)
            end

        end

        function saveData(obj,step,cParams)
            obj.data.displacement.function{step} = cParams.uFun;
            obj.data.displacement.value(step)    = cParams.uVal;
            obj.data.reaction                    = cParams.fVal;
            obj.data.damage.function{step}       = cParams.dmgFun;
            obj.data.damage.maxValue(step)       = max(cParams.dmgFun.fValues);
            obj.data.qMaxValue                   = cParams.qMax;
            obj.data.rMaxValue                   = cParams.rMax;
            obj.data.energy(step)                = cParams.energy;
            obj.data.iter(step)                  = cParams.numIter;
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            s.shallDisplay = cParams.shallDisplay;
            s.maxNColumns = 3;
            s.funs = [{cParams.fun}];
            s.barLims = [{[0;1]}];
            s.titles = [{'Force-displacement'},{'Damage-displacement'},{'q-r'},{'Cost'},{'Damage'},{'Iterations'}];
            s.chartTypes = [{'plot'},{'plot'},{'plot'},{'plot'},{'surf'},{'bar'}];
            obj.monitor = Monitoring(s);

            obj.data = [];
            obj.printInfo = cParams.shallPrintInfo;
        end

    end
    
end