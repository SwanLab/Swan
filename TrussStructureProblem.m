classdef TrussStructureProblem < handle
    
    properties (Access = public)
        meshFileName = 'ISCSOMesh'
        barLengthFileName = 'barLength.mat'
        result
    end

    properties (Access = private)
        
    end

    methods (Access = public)
        
        function obj = TrussStructureProblem()
            obj.compute();
        end

    end
       
    methods (Access = private)
        
        function compute(obj)
            cParams   = obj.getDesignVarParams();
            designVar = DesignVariableTruss(cParams);
            nBars     = size(connec,1);
            x0        = [ones(nBars,1); ones(nBars,1)];
            designVar.init(x0);
            interp             = BarSectionInterpolation(designVar);
            costData.designVar = designVar;
            costData.interp    = interp;
            costData.barLength = obj.barLengthFileName;
            constrData.phyProb = obj.createFEMProblem(interp, designVar);
            constrData.interp  = interp;
            s.designVar        = designVar;
            s.cost             = TrussStructureCost(costData);
            s.constraint       = TrussStructureConstraint(constrData);
            s.outputFunction.type       = "Truss Structure";
            s.outputFunction.monitoring = MonitoringManager(s);
            s.optimizerNames.primal     = 'PROJECTED GRADIENT';
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
            prob = StructuralTrussProblem(s);
        end

    end
    
end