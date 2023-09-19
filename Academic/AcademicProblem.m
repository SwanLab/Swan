classdef AcademicProblem < handle

    properties (Access = public)
        result
    end
    
    properties (Access = private)
        cost
        constraint
        optimizer
        filename
    end

    methods (Access = public)

        function obj = AcademicProblem(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.filename = cParams.filename;
        end

        function createOptimizer(obj)
            run(obj.filename);
            p.x0                         = x0;
            d                            = DesignVariableAcademic(p);
            j.dV                         = d;
            c.dV                         = d;
            s.designVar                  = d;
            s.cost                       = AcademicCost(j);
            s.constraint                 = AcademicConstraint(c);
            s.constraint.nSF             = nConstr;
            s.nConstraints               = nConstr;
            s.dualVariable               = DualVariable(s);
            s.outputFunction.type        = "Academic";
            s.outputFunction.iterDisplay = "iter";
            s.shallPrint = shallPrint;
            s.outputFunction.monitoring  = MonitoringManager(s);
            s.optimizerNames.primal     = 'PROJECTED GRADIENT';
            opt = Optimizer.create(s);
            opt.solveProblem();
            obj.result = d.value;
        end

    end

end