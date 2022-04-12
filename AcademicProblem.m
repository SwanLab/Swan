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
            s.constraintCase             = [];
            s.designVar                  = d;
            s.dualVariable               = [];
            s.maxIter                    = [];
            s.incrementalScheme          = [];
            s.targetParameters           = [];
            s.cost                       = AcademicCost(j);
            s.constraint                 = AcademicConstraint(c);
            s.outputFunction.type        = "Academic";
            s.outputFunction.iterDisplay = "iter";
            opt = Optimizer.create(s);
            opt.solveProblem();
            d.value
        end

    end

end