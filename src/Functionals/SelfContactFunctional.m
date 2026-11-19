classdef SelfContactFunctional < BoundaryContactFunctional

    properties (Access = protected)
        penalty
        depth
    end

    methods (Access = public)

        function obj = SelfContactFunctional(cParams)
            obj@BoundaryContactFunctional(cParams);
            obj.penalty = cParams.penalty;
            obj.depth   = cParams.depth;
        end

    end

    methods (Access = protected)

        function c = computeContacts(obj,uFun)
            x  = obj.computeCurrentCoordinates(uFun);
            [node,edge] = obj.computeCandidatePairs(x,obj.depth);
            a  = obj.edges(edge,1);
            b  = obj.edges(edge,2);
            [xi,g]       = obj.computeProjection(x,node,a,b);
            [xiRef,gRef] = obj.computeProjection(obj.mesh.coord,node,a,b);
            isBehindInReference = gRef <= 0 & xiRef > -1 & xiRef < 2;
            active = node ~= a & node ~= b & ~isBehindInReference & ...
                     xi > -obj.windowOverhang() & xi < 1+obj.windowOverhang() & ...
                     obj.isInContactRange(node,g);
            obj.storeActivePairs(node(active));
            c = obj.packPairs(x,node(active),a(active),b(active));
        end

        function active = isInContactRange(obj,~,g)
            active = g < 0 & g > -obj.depth;
        end

        function storeActivePairs(~,~)
        end

        function e = computePairEnergy(obj,X)
            r  = X(:,1:2) - X(:,3:4);
            d  = X(:,5:6) - X(:,3:4);
            xi = sum(r.*d,2)./sum(d.^2,2);
            g  = (r(:,1).*d(:,2) - r(:,2).*d(:,1))./sqrt(sum(d.^2,2));
            e  = obj.penalty/2*obj.computeWindow(xi).*g.^2;
        end

        function v = computePairGradient(obj,X)
            r  = X(:,1:2) - X(:,3:4);
            d  = X(:,5:6) - X(:,3:4);
            L2 = sum(d.^2,2);
            xi = sum(r.*d,2)./L2;
            n  = [d(:,2),-d(:,1)]./sqrt(L2);
            g  = sum(r.*n,2);
            [w,dw] = obj.computeWindow(xi);
            dGap = [n,-(1-xi).*n,-xi.*n];
            dXi  = [d,-r+(2*xi-1).*d,r-2*xi.*d]./L2;
            v = obj.penalty*(w.*g.*dGap + 0.5*dw.*g.^2.*dXi);
        end

    end

    methods (Access = private)

        function [xi,g] = computeProjection(~,x,node,a,b)
            d  = x(b,:) - x(a,:);
            r  = x(node,:) - x(a,:);
            xi = sum(r.*d,2)./sum(d.^2,2);
            n  = [d(:,2),-d(:,1)]./sqrt(sum(d.^2,2));
            g  = sum(r.*n,2);
        end

    end

    methods (Static, Access = protected)

        function [w,dw] = computeWindow(xi)
            e = SelfContactFunctional.windowOverhang();
            s = min(max(min(xi+e,1+e-xi)/e,0),1);
            w = 3*s.^2 - 2*s.^3;
            ds = (xi < 0).*(xi > -e)/e - (xi > 1).*(xi < 1+e)/e;
            dw = 6*s.*(1-s).*ds;
        end

        function e = windowOverhang()
            e = 0.2;
        end

    end

end
