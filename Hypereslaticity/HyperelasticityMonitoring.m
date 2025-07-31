classdef HyperelasticityMonitoring < handle

    properties (GetAccess = public, SetAccess = private)
        data
    end

    properties (Access = private)
        monitor
        printInfo
        printFile
        fileNameOut
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

        function printOutput(obj,step,u,r)
            if obj.printFile
                fun = {u,r};
                funNames = {'Displacement', 'Reactions'};
                a.mesh     = u.mesh;
                a.filename = ['SIM_',obj.fileNameOut,'_',int2str(step)];
                a.fun      = fun;
                a.funNames = funNames;
                a.type     = 'Paraview';
                pst = FunctionPrinter.create(a);
                pst.print();
            end
        end

        function saveData(obj,step,cParams)
            obj.data.reaction.value(step)        = cParams.r;
            obj.data.reaction.function{step}     = cParams.rFun;
            obj.data.displacement.value(step)    = cParams.u;
            obj.data.displacement.function{step} = cParams.uFun;
            obj.data.energy.linear(step)         = cParams.energy;
            obj.data.iter(step)                  = cParams.numIter;
        end
        
    end

    methods (Access = private)
        
        function init(obj,cParams)
            s.shallDisplay = cParams.shallDisplay;
            s.maxNColumns = 2;
            s.titles = [{'Energy'},{'Cost'},{'Force-displacement'},{'Iterations'}];
            s.chartTypes = [{'plot'},{'plot'},{'plot'},{'bar'}];
            obj.monitor = Monitoring(s);

            obj.data = [];
            obj.printInfo = cParams.shallPrintInfo;
            obj.printFile = cParams.shallPrintFile;
            obj.fileNameOut = cParams.fileNameOut;
        end

    end
    
end