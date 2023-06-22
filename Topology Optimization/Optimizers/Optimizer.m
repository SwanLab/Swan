classdef Optimizer < handle
    
    properties (Access = protected)
        designVariable
        dualVariable
        cost
        constraint
        outputFunction
        maxIter
        nIter = 0
        targetParameters
        dualUpdater
        primalUpdater
        constraintCase
        postProcess
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        type
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
        
    end
    
    methods (Access = protected)
        
        function initOptimizer(obj,cParams)
            obj.nIter             = 0;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.maxIter           = cParams.maxIter;
            obj.targetParameters  = cParams.targetParameters;
            obj.constraintCase    = cParams.constraintCase;
            obj.outputFunction    = cParams.outputFunction.monitoring;
            obj.createPostProcess(cParams.postProcessSettings);
        end

        function createPrimalUpdater(obj,cParams)
            f                 = PrimalUpdaterFactory();
            obj.primalUpdater = f.create(cParams);
        end

        function createDualUpdater(obj,cParams)
            f               = DualUpdaterFactory();
            obj.dualUpdater = f.create(cParams);      
        end
        
        function isAcceptable = checkConstraint(obj)
            for i = 1:length(obj.constraint.value)
                switch obj.constraintCase{i}
                    case {'EQUALITY'}
                        fine = obj.checkEqualityConstraint(i);
                    case {'INEQUALITY'}
                        fine = obj.checkInequalityConstraint(i);
                end
                if fine
                    isAcceptable = true;
                else
                    isAcceptable = false;
                end
            end
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
%                 obj.simulationPrinter.appendStep(file);
            end
        end

    end

    methods (Access = private)

        function createPostProcess(obj,cParams)
            if cParams.shallPrint
                d = obj.createPostProcessDataBase(cParams);
                d.printMode = cParams.printMode;
                d.nDesignVariables = obj.designVariable.nVariables;
                obj.postProcess = Postprocess('TopOptProblem',d);
                s.filename = [obj.outFolder,'/',obj.outFilename, '_simulation'];
%                 obj.simulationPrinter = SimulationPrinter(s);
            end
        end

        function d = createPostProcessDataBase(obj,cParams)
            path = pwd;
            obj.outFolder   = fullfile(path,'Output',cParams.femFileName);
            obj.outFilename = cParams.femFileName;
            d.mesh    = obj.designVariable.mesh;
            d.outFileName = cParams.femFileName;
            d.ptype   = cParams.ptype;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.create();
            d.ndim       = obj.designVariable.mesh.ndim;
            d.pdim       = cParams.pdim;
            d.optimizer  = obj.type;
            d.cost       = obj.cost;
            d.constraint = obj.constraint;
            d.designVar  = obj.designVariable.type;
        end

        function c = checkInequalityConstraint(obj,i)
            g = obj.constraint.value(i);
            c = g <= 0;
        end

        function c = checkEqualityConstraint(obj,i)
            g = obj.constraint.value(i);
            c = abs(g) < obj.targetParameters.constr_tol;
        end

    end

end