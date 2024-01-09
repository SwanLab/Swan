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
            obj.compute();
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
            s.outputFunction.iterDisplay = "off";
            s.outputFunction.monitoring  = MonitoringManager(s);
            s.optimizerNames.primal     = 'PROJECTED GRADIENT';
            opt = Optimizer.create(s);
            opt.solveProblem();
            obj.result = d.value;
            close all;
        end

    end

end