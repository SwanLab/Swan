classdef AugmentedLagrangianContactFunctional < SelfContactFunctional

    properties (Access = private)
        multipliers
        pairMultipliers
        penetrationTolerance
    end

    methods (Access = public)

        function obj = AugmentedLagrangianContactFunctional(cParams)
            obj@SelfContactFunctional(cParams);
            obj.multipliers          = zeros(cParams.mesh.nnodes,1);
            obj.penetrationTolerance = cParams.penetrationTolerance;
        end

        function itHas = updateParameters(obj,uFun)
            c = obj.computeContacts(uFun);
            [~,g,node] = obj.computePairGaps(c.coords,c.dofs);
            newMultipliers = zeros(size(obj.multipliers));
            if ~isempty(node)
                updated = max(0,obj.pairMultipliers - obj.penalty*g);
                newMultipliers = accumarray(node,updated,size(obj.multipliers),@max);
            end
            maxPenetration = max([0;-g]);
            obj.multipliers = newMultipliers;
            itHas = maxPenetration < obj.penetrationTolerance;
        end

    end

    methods (Access = protected)

        function active = isInContactRange(obj,node,g)
            active = g < obj.multipliers(node)/obj.penalty & g > -obj.depth;
        end

        function storeActivePairs(obj,node)
            obj.pairMultipliers = obj.multipliers(node);
        end

        function e = computePairEnergy(obj,X)
            [xi,g] = obj.computePairGaps(X);
            lambda = obj.pairMultipliers;
            e = obj.computeWindow(xi).*(max(0,lambda - obj.penalty*g).^2 - lambda.^2)/(2*obj.penalty);
        end

        function v = computePairGradient(obj,X)
            r  = X(:,1:2) - X(:,3:4);
            d  = X(:,5:6) - X(:,3:4);
            L2 = sum(d.^2,2);
            xi = sum(r.*d,2)./L2;
            n  = [d(:,2),-d(:,1)]./sqrt(L2);
            g  = sum(r.*n,2);
            lambda = obj.pairMultipliers;
            [w,dw] = obj.computeWindow(xi);
            pressure = max(0,lambda - obj.penalty*g);
            dGap = [n,-(1-xi).*n,-xi.*n];
            dXi  = [d,-r+(2*xi-1).*d,r-2*xi.*d]./L2;
            v = -w.*pressure.*dGap + dw.*(pressure.^2 - lambda.^2)/(2*obj.penalty).*dXi;
        end

    end

    methods (Access = private)

        function [xi,g,node] = computePairGaps(~,X,dofs)
            r  = X(:,1:2) - X(:,3:4);
            d  = X(:,5:6) - X(:,3:4);
            xi = sum(r.*d,2)./sum(d.^2,2);
            g  = (r(:,1).*d(:,2) - r(:,2).*d(:,1))./sqrt(sum(d.^2,2));
            node = [];
            if nargin > 2
                node = (dofs(:,2))/2;
            end
        end

    end

end
