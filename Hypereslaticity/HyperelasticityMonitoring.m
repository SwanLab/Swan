classdef HyperelasticityMonitoring < handle

    properties (GetAccess = public, SetAccess = private)
        data
    end

    properties (Access = private)
        monitor
        print
    end

    methods (Access = public)

        function obj = HyperelasticityMonitoring(cParams)
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
            if obj.print == true
                X = sprintf('%s:%d / cost: %.8e  (diff:%.8e) \n',name,iter,cost,err);
                fprintf(X);
            end
        end

        function printStep(obj,step,maxSteps)
            if obj.print == true
                fprintf('\n ********* STEP %i/%i *********  \n',step,maxSteps)
            end

        end

        function saveData(obj,step,cParams)
            obj.data.reaction.value(step)        = cParams.r;
            obj.data.reaction.function{step}     = cParams.rFun;
            obj.data.displacement.value(step)    = cParams.u;
            obj.data.displacement.function{step} = cParams.uFun;

            obj.data.energy.linear(step)     = cParams.energy(1);
            obj.data.energy.neohook(step)    = cParams.energy(2);

            obj.data.iter(step)    = cParams.numIter;
        end
        
    end

    methods (Access = private)
        
        function init(obj,cParams)
            s.shallDisplay = cParams.shallDisplay;
            s.maxNColumns = 3;
            s.titles = [{'Energy'},{'Force-displacement'},{'Iterations'}];
            s.chartTypes = [{'plot'},{'plot'},{'bar'}];
            obj.monitor = Monitoring(s);

            obj.data = [];
            obj.print = cParams.shallPrint;
        end

    end
    
end