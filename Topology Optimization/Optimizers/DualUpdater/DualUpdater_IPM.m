classdef DualUpdater_IPM < handle

    properties (Access = private)
        dualVariable
        constraint
        nConstr
        constraintCase
        cost
        alpha
    end

    methods (Access = public)
        
        function obj = DualUpdater_IPM(cParams)
            obj.init(cParams);
        end

        function compute(obj,lz,uz)
            c   = obj.constraint.gradient';
            g   = obj.cost.gradient;
            l   = (c*c')\(c*(lz' - uz' - g'));
            obj.dualVariable.value = l';
        end

        function updateAlpha(obj,a)
            obj.alpha = a;
        end

        function obj = update(obj,g)
            tau                    = obj.alpha;
            obj.dualVariable.value = obj.dualVariable.value + tau * g';
        end

        function zLB = updateLowerBound(obj,bounds)
            tau  = obj.alpha;
            zLB  = bounds.zLB;
            dzLB = bounds.dzLB;
            zLB  = zLB + tau*dzLB';
        end

        function zUB = updateUpperBound(obj,bounds)
            tau  = obj.alpha;
            zUB  = bounds.zUB;
            dzUB = bounds.dzUB;
            zUB  = zUB + tau*dzUB';
        end
    end


    methods (Access = private)

        function init(obj,cParams)
            obj.dualVariable   = cParams.dualVariable;
            obj.constraint     = cParams.constraint;
            obj.constraintCase = cParams.constraintCase;
            obj.nConstr        = cParams.constraint.nSF;
            obj.cost           = cParams.cost;
        end

    end

    methods (Static, Access = public)
        function zLB = computeLowerBound(cParams)
            mu  = cParams.mu;
            x   = cParams.designVariable.fun.fValues';
            xLB = cParams.bounds.xLB;
            s   = cParams.slack';
            sLB = cParams.bounds.sLB;
            zLB = mu./(x-xLB);
            zLB = [zLB,mu./(s-sLB)];
        end

        function zUB = computeUpperBound(cParams)
            mu  = cParams.mu;
            x   = cParams.designVariable.fun.fValues';
            xUB = cParams.bounds.xUB;
            s   = cParams.slack';
            sUB = cParams.bounds.sUB;
            zUB = mu./(xUB-x);
            zUB = [zUB,mu./(sUB-s)];
        end
    end

end