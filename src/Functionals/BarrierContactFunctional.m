classdef BarrierContactFunctional < BoundaryContactFunctional

    properties (Access = private)
        stiffness
        initialStiffness
        stiffnessReduction
        gapTolerance
        activationDistance
    end

    methods (Access = public)

        function obj = BarrierContactFunctional(cParams)
            obj@BoundaryContactFunctional(cParams);
            obj.initialStiffness   = cParams.stiffness;
            obj.stiffness          = cParams.stiffness;
            obj.activationDistance = cParams.activationDistance;
            obj.initContinuation(cParams);
        end

        function resetParameters(obj)
            obj.stiffness = obj.initialStiffness;
        end

        function itHas = updateParameters(obj,uFun)
            itHas = obj.isGapClosed(uFun);
            if ~itHas
                obj.stiffness = obj.stiffness/obj.stiffnessReduction;
            end
        end

        function cost = computeCost(obj,uFun)
            cost = computeCost@BoundaryContactFunctional(obj,uFun);
            if obj.countOverlaps(uFun) > 0
                cost = Inf;
            end
        end

    end

    methods (Access = protected)

        function c = computeContacts(obj,uFun)
            x    = obj.computeCurrentCoordinates(uFun);
            dHat = obj.activationDistance;
            [node,edge] = obj.computeCandidatePairs(x,dHat);
            a = obj.edges(edge,1);
            b = obj.edges(edge,2);
            d    = obj.computeDistance(obj.gatherPairCoordinates(x,node,a,b));
            dRef = obj.computeDistance(obj.gatherPairCoordinates(obj.mesh.coord,node,a,b));
            active = node ~= a & node ~= b & dRef > 2*dHat & d < dHat;
            c = obj.packPairs(x,node(active),a(active),b(active));
        end

        function e = computePairEnergy(obj,X)
            dHat = obj.activationDistance;
            d = obj.computeDistance(X);
            e = -obj.stiffness*(d-dHat).^2.*log(d/dHat);
        end

        function v = computePairGradient(obj,X)
            dHat = obj.activationDistance;
            [d,dd] = obj.computeDistance(X);
            db = -2*(d-dHat).*log(d/dHat) - (d-dHat).^2./d;
            v  = obj.stiffness*db.*dd;
        end

    end

    methods (Access = private)

        function itIs = isGapClosed(obj,uFun)
            c    = obj.computeContacts(uFun);
            itIs = isempty(c.coords) || min(obj.computeDistance(c.coords)) < obj.gapTolerance;
        end

        function initContinuation(obj,cParams)
            obj.gapTolerance       = Inf;
            obj.stiffnessReduction = 1;
            if isfield(cParams,'gapTolerance')
                obj.gapTolerance       = cParams.gapTolerance;
                obj.stiffnessReduction = cParams.stiffnessReduction;
            end
        end

    end

    methods (Static, Access = private)

        function [dist,dDist] = computeDistance(X)
            p  = X(:,1:2);
            a  = X(:,3:4);
            r  = p - a;
            t  = X(:,5:6) - a;
            xi = min(max(sum(r.*t,2)./sum(t.^2,2),0),1);
            v  = r - xi.*t;
            dist  = sqrt(sum(v.^2,2));
            n     = v./dist;
            dDist = [n,-(1-xi).*n,-xi.*n];
        end

    end

end
