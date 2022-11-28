classdef TrussStructureProblem < handle
    
    properties (Access = public)
        fileName = 'DataDoc'
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
            designVar          = DesignVariableTruss(cParams);
            designVar.init(x0);
            interp             = BarSectionInterpolation(d);
            costData.designVar = designVar;
            costData.interp    = interp;
            constrData.phyProb = TrussStructureFEM(); % AQUI EL QUE NECESSITIS TON QUE ENTRI JA DIRÀS
            s.designVar        = designVar;
            s.cost             = TrussStructureCost(costData);
            s.constraint       = TrussStructureConstraint(constrData);
        

        end

    end
    
end