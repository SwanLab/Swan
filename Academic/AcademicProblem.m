classdef AcademicProblem

    properties (Access = private)
        cost
        constraint
        optimizer
    end

    methods (Access = public)

        function obj = AcademicProblem()
            obj.init();
            obj.compute();
        end

        function compute(obj)
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            
        end

        function createOptimizer(obj)
            run("AcademicTest1.m");
            d                            = DesignVariableAcademic();
            d.init(x0);
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
            s.outputFunction.monitoring  = MonitoringManager(s);
            s.optimizerNames.primal     = 'PROJECTED GRADIENT';
            opt = Optimizer.create(s);
            opt.solveProblem();
            d.value
        end

    end

end