classdef MonitoringOptimizationProblemStandard < handle

    properties (Access = public)
        nPlots
    end

    properties (Access = private)
        cost
        constraint
        designVariable
        dualVariable
        isConstrained
    end

    methods (Access = public)
        function obj = MonitoringOptimizationProblemStandard(cParams)
            obj.init(cParams);
        end

        function m = create(obj,m)
            m = obj.createCostMonitoring(m);
            m = obj.createConstraintsMonitoring(m);
            m = obj.createDesignVariableNormMonitoring(m);
            m = obj.createLagrangeMultipliersMonitoring(m);
        end

        function m = plot(obj,m,it)
            m = obj.plotCost(m,it);
            m = obj.plotConstraints(m,it);
            m = obj.plotDesignVariableNorm(m,it);
            m = obj.plotLagrangeMultipliers(m,it);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.designVariable = cParams.designVariable;
            obj.isConstrained  = cParams.isConstrained;
            switch cParams.isConstrained
                case true
                    nConstr          = cParams.constraint.obtainNumberFields();
                    obj.constraint   = cParams.constraint;
                    obj.dualVariable = cParams.dualVariable;
                case false
                    nConstr  = 0;
            end
            nF         = obj.cost.obtainNumberFields();
            obj.nPlots = 2+nF+2*nConstr;
        end

        function m = createCostMonitoring(obj,m)
            titlesF = obj.cost.getTitleFields();
            m       = MonitoringVariable.create(m,1,'Cost',obj.cost);
            for i = 1:length(titlesF)
                m = MonitoringVariable.create(m,1+i,titlesF{i},obj.cost);
            end
        end

        function m = createConstraintsMonitoring(obj,m)
            if obj.isConstrained
                titlesConst = obj.constraint.getTitleFields();
                n           = length(m.figures);
                for i = 1:length(titlesConst)
                    m = MonitoringVariable.create(m,n+i,titlesConst{i},obj.constraint);
                end
            end
        end

        function m = createDesignVariableNormMonitoring(obj,m)
            n = length(m.figures);
            m = MonitoringVariable.create(m,n+1,'Norm L2 x',obj.designVariable);
        end

        function m = createLagrangeMultipliersMonitoring(obj,m)
            if obj.isConstrained
                titlesConst = obj.constraint.getTitleFields();
                n           = length(m.figures);
                for i = 1:length(titlesConst)
                    titleLambda = ['\lambda_{',titlesConst{i},'}'];
                    m = MonitoringVariable.create(m,n+1,titleLambda,obj.dualVariable);
                end
            end
        end

        function m = plotCost(obj,m,it)
            nF = obj.cost.obtainNumberFields();
            m.figures{1}.updateParams(it,m.data{1}.value);
            m.figures{1}.refresh();
            for i = 1:nF
                m.figures{i+1}.updateParams(it,m.data{i+1}.getFields(i));
                m.figures{i+1}.refresh();
            end
        end

        function m = plotConstraints(obj,m,it)
            nF      = obj.cost.obtainNumberFields();
            n       = nF+1;
            nConstr = obj.constraint.obtainNumberFields();
            for i = 1:nConstr
                m.figures{i+n}.updateParams(it,m.data{i+n}.value(i,1));
                m.figures{i+n}.refresh();
            end
        end

        function m = plotDesignVariableNorm(obj,m,it)
            nF      = obj.cost.obtainNumberFields();
            nConstr = obj.constraint.obtainNumberFields();
            n       = nF+nConstr+1;
            m.figures{1+n}.updateParams(it,m.data{1+n}.computeL2normIncrement());
            m.figures{1+n}.refresh();
        end

        function m = plotLagrangeMultipliers(obj,m,it)
            nF      = obj.cost.obtainNumberFields();
            nConstr = obj.constraint.obtainNumberFields();
            n       = nF+nConstr+2;
            for i = 1:nConstr
                m.figures{i+n}.updateParams(it,m.data{i+n}.value(i,1));
                m.figures{i+n}.refresh();
            end
        end
    end
end