classdef TrussStructureProblem < handle
    
    properties (Access = public)
        meshFileName = 'ISCSOMesh'
        problemFile  = 'ISCSOProblem'
        result
    end

    properties (Access = private)
        barsLength
        nBars
    end

    methods (Access = public)
        
        function obj = TrussStructureProblem()
            obj.createAreaAndInertia();
            obj.compute();
        end

    end
       
    methods (Access = private)
        
        function compute(obj)
            cParams   = obj.getDesignVarParams();
            run(obj.problemFile);
            designVar = DesignVariableTruss(cParams);
            nBars = 6;
            x0 = [ones(nBars,1); ones(nBars,1)];
            designVar.init(x0);
            interp              = BarSectionInterpolation(designVar);
            costData.designVariable = designVar;
            costData.interp     = interp;
            costData.barsLength = ones(6,1);
            constrData.phyProb  = obj.createFEMProblem(interp, designVar);
            constrData.interp   = interp;
            constrData.nConstraints = nConstr;
            constrData.designVariable = designVar;
            s.designVar         = designVar;
            s.nConstraints      = nConstr;
            s.cost              = TrussStructureCost(costData);
            s.constraint        = TrussStructureConstraint(constrData);
            s.constraint.nSF    = nConstr;
            s.dualVariable      = DualVariable(s);
            s.type              = 'NullSpace';
            s.outputFunction.type       = "Academic";
            s.incrementalScheme.nSteps  = 1;
            s.outputFunction.monitoring = MonitoringManager(s);
            s.optimizerNames.primal     = 'PROJECTED GRADIENT';
            s.maxIter                   = 1000;
            s.constraintCase    = {'INEQUALITY'};
            s.targetParameters  = 0;
            s.postProcessSettings.shallPrint = false;
            opt = Optimizer.create(s);
            opt.solveProblem();
            obj.result = designVar.value;
        end

        function s = getDesignVarParams(obj)
            s.ubR = 0;
            s.lbR = 0;
            s.ubT = 0;
            s.lbT = 0;
        end

        function prob = createFEMProblem(obj, interp, dVar)
            run(obj.meshFileName); % AQUÍ HA D'HAVER-HI TOTA LA INFO
            s.coord     = coord;
            s.connec    = connec;
            s.material  = ISCSOMaterial(s);
            s.neumann   = neumann;
            s.dirichlet = dirichlet;
            s.interp    = interp;
            s.designVar = dVar;
            obj.nBars   = size(connec,1);
            obj.barsLength = ones(6,1);% S'HA DE CALCULAR
            prob = StructuralTrussProblem(s);
        end

        function createAreaAndInertia()
            Si = DataSections();
            A = Si(:,1);
            I = Si(:,2);
        end

    end
    
end