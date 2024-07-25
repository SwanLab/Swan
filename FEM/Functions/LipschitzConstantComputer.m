classdef LipschitzConstantComputer < handle

    properties (Access = private)
        cost
        designVariable
        isLagrangian
        constraint
    end

    methods (Access = public)
        function obj = LipschitzConstantComputer(cParams)
            obj.init(cParams);
        end

        function L = compute(obj,sampleX)
            nSample = size(sampleX,2);
            L       = 0;
            for i = 1:nSample
                xi     = sampleX(:,i);
                x      = obj.designVariable;
                x.update(xi);
                [~,dJ] = obj.cost.computeFunctionAndGradient(x);
                if obj.isLagrangian
                    [~,dg]    = obj.constraint{1}.computeFunctionAndGradient(x);
                    lEstimate = sum(abs(dJ.fValues))./sum(abs(dg.fValues),1);
                    dM        = dJ.fValues + dg.fValues*lEstimate';
                else
                    dM = dJ;
                    dg = 0;
                end
                L = max(L,norm(dJ.fValues));
                L = max(L,norm(dg.fValues));
                L = max(L,norm(dM));
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.designVariable = cParams.designVariable.copy();
            if isfield(cParams,'constraint')
                obj.isLagrangian      = true;
                obj.constraint        = cParams.constraint;
            else
                obj.isLagrangian = false;
            end
        end
    end
end