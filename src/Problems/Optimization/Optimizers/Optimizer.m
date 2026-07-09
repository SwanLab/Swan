classdef Optimizer < handle

    properties (Access = protected)
        designVariable
        dualVariable
        cost
        constraint
        monitoring
        maxIter
        nIter = 0
        tolerance
        dualUpdater
        primalUpdater
        constraintCase
        postProcess
    end

    properties (Access = public)
        simulationPrinter
    end

    properties (Access = private)
        outFilename % !!!
        outFolder
    end


    methods (Access = public, Static)

        function obj = create(cParams)
            f   = OptimizerFactory();
            obj = f.create(cParams);
        end

        function isAcceptable = checkConstraint(value,cases,tol)
            for i = 1:length(value)
                switch cases{i}
                    case {'EQUALITY'}
                        fine = Optimizer.checkEqualityConstraint(value,i,tol);
                    case {'INEQUALITY'}
                        fine = Optimizer.checkInequalityConstraint(value,i);
                end
                if fine
                    isAcceptable = true;
                else
                    isAcceptable = false;
                    break;
                end
            end
        end

    end

    methods (Access = protected)

        function initOptimizer(obj,cParams)
            obj.nIter          = 0;
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.designVariable = cParams.designVariable;
            obj.dualVariable   = cParams.dualVariable;
            obj.maxIter        = cParams.maxIter;
            obj.tolerance      = cParams.tolerance;
            obj.constraintCase = cParams.constraintCase;
        end

        function createPrimalUpdater(obj,cParams)
            f                 = PrimalUpdaterFactory();
            obj.primalUpdater = f.create(cParams);
        end

        function createDualUpdater(obj,cParams)
            cParams.type    = obj.type;
            f               = DualUpdaterFactory();
            obj.dualUpdater = f.create(cParams);
        end

        function printOptimizerVariable(obj)
            if ~isempty(obj.postProcess)
                d.fields  = obj.designVariable.getVariablesToPlot();
                d.cost = obj.cost;
                d.constraint = obj.constraint;
                %                 obj.postProcess.print(obj.nIter,d);
                [desFun, desName] = obj.designVariable.getFunsToPlot();
                fun  = desFun;
                name = desName;
                for iShp = 1:numel(obj.cost.shapeFunctions)
                    [shpFun, shpName] = obj.cost.shapeFunctions{iShp}.getFunsToPlot();
                    fun  = [fun, shpFun];
                    name = [name, shpName];
                end
                file = [obj.outFolder,'/',obj.outFilename, '_', num2str(obj.nIter)];

                zz.mesh     = obj.designVariable.mesh;
                zz.filename = file;
                zz.fun      = fun;
                zz.funNames = name;
                pp = FunctionPrinter_Paraview(zz);
                pp.print();
                obj.simulationPrinter.appendStep(file);
            end
            if ismethod(obj.designVariable,'plot')
                obj.designVariable.plot();
            end
        end
    end

    methods (Static, Access = private)
        function c = checkInequalityConstraint(value,i)
            g = value(i);
            c = g <= 0;
        end

        function c = checkEqualityConstraint(value,i,tol)
            g = value(i);
            c = abs(g) < tol;
        end
    end
end