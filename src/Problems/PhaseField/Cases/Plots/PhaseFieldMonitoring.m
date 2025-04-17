classdef PhaseFieldMonitoring < handle


    properties (GetAccess = public, SetAccess = private)
        data
    end

    properties (Access = private)
        monitor
        print
    end

    methods (Access = public)

        function obj = PhaseFieldMonitoring(cParams)
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
            obj.data.force(step)              = cParams.force;
            obj.data.displacement.value(step) = cParams.bcVal;
            obj.data.displacement.field       = cParams.u;
            obj.data.damage.field             = cParams.phi;
            obj.data.damage.maxValue(step)    = max(cParams.phi.fValues);

            obj.data.energy.intE(step)     = cParams.energy(1);
            obj.data.energy.localDis(step) = cParams.energy(2);
            obj.data.energy.regDis(step)   = cParams.energy(3);
            obj.data.energy.extWork(step)  = cParams.energy(4);

            obj.data.iter.u(step)    = cParams.numIterU;
            obj.data.iter.phi(step)  = cParams.numIterP;
            obj.data.iter.stag(step) = cParams.numIterStag;

            obj.data.cost = cParams.cost;
            obj.data.tau  = cParams.tauArray;
        end
        
    end

    methods (Access = private)
        
        function init(obj,cParams)
            s.shallDisplay = cParams.shallDisplay;
            s.funs = [{cParams.fun}];
            s.barLims = [{[0;1]}];
            switch cParams.type
                case 'full'
                    s.maxNColumns = 3;
                    s.titles = [{'Force-displacement'},{'Damage-displacement'},{'Damage'},{'Iter Staggered'},{'Cost'},{'Total Energy'},{'Line-search'}]; % ,{'Energy'}
                    s.chartTypes = [{'plot'},{'plot'},{'surf'},{'plot'},{'plot'},{'plot'},{'plot'}];
                    obj.monitor = Monitoring(s);
                case 'reduced'
                    s.maxNColumns = 2;
                    s.titles = [{'Force-displacement'},{'Damage-displacement'},{'Damage'}];
                    s.chartTypes = [{'plot'},{'plot'},{'surf'}];
                    obj.monitor = Monitoring(s);
            end
            obj.data = [];
            obj.print = cParams.shallPrint;
        end

    end
    
end