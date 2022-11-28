classdef TrussStructureProblem < handle
    
    properties (Access = public)
        fileName = 'ISCSOMesh'
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
            run(obj.fileName); % AQUÍ HA D'HAVER-HI TOTA LA INFO
            cParams = obj.getDesignVarParams();
            designVar          = DesignVariableTruss(cParams);
            nBars = size(connec,1);
            x0 = [ones(nBars,1); ones(nBars,1)];
            designVar.init(x0);
            interp             = BarSectionInterpolation(designVar);
            costData.designVar = designVar;
            costData.interp    = interp;
            constrData.phyProb = obj.createFEMProblem(interp, designVar);
            s.designVar        = designVar;
            s.cost             = TrussStructureCost(costData);
            s.constraint       = TrussStructureConstraint(constrData);
        

        end

        function s = getDesignVarParams(obj)
            s.ubR = 0;
            s.lbR = 0;
            s.ubT = 0;
            s.lbT = 0;
        end

        function prob = createFEMProblem(obj, interp, dVar)
            run(obj.fileName); % AQUÍ HA D'HAVER-HI TOTA LA INFO
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