classdef BoundaryContactFunctional < handle

    properties (Access = protected)
        mesh
        edges
        boundaryNodes
    end

    methods (Access = public)

        function obj = BoundaryContactFunctional(cParams)
            obj.mesh = cParams.mesh;
            obj.computeOrientedBoundaryEdges();
        end

        function cost = computeCost(obj,uFun)
            c = obj.computeContacts(uFun);
            cost = sum(obj.computePairEnergy(c.coords));
        end

        function grad = computeGradient(obj,uFun)
            c = obj.computeContacts(uFun);
            v = obj.computePairGradient(c.coords);
            grad = accumarray(c.dofs(:),v(:),[obj.nDofs() 1]);
        end

        function hess = computeHessian(obj,uFun)
            c = obj.computeContacts(uFun);
            k = obj.computePairHessian(c.coords);
            i = repmat(c.dofs,[1 1 6]);
            j = permute(i,[1 3 2]);
            hess = sparse(i(:),j(:),k(:),obj.nDofs(),obj.nDofs());
        end

        function resetParameters(~)
        end

        function itHas = updateParameters(~,~)
            itHas = true;
        end

        function n = countOverlaps(obj,uFun)
            x = obj.computeCurrentCoordinates(uFun);
            p = x(obj.edges(:,1),:);
            q = x(obj.edges(:,2),:);
            h = max(sqrt(sum((q-p).^2,2)));
            [i,j] = obj.findNeighbours((p+q)/2,(p+q)/2,h);
            ei = obj.edges(i,:);
            ej = obj.edges(j,:);
            shareNode = ej(:,1) == ei(:,1) | ej(:,1) == ei(:,2) | ej(:,2) == ei(:,1) | ej(:,2) == ei(:,2);
            keep = j > i & ~shareNode;
            i = i(keep);
            j = j(keep);
            n = sum(obj.doSegmentsCross(p(i,:),q(i,:),p(j,:),q(j,:)));
        end

    end

    methods (Abstract, Access = protected)
        c = computeContacts(obj,uFun)
        e = computePairEnergy(obj,X)
        v = computePairGradient(obj,X)
    end

    methods (Access = protected)

        function x = computeCurrentCoordinates(obj,uFun)
            x = obj.mesh.coord + uFun.fValues(1:obj.mesh.nnodes,:);
        end

        function c = packPairs(obj,x,node,a,b)
            nodes    = [node,a,b];
            c.coords = obj.gatherPairCoordinates(x,node,a,b);
            dofs     = [2*nodes-1,2*nodes];
            c.dofs   = dofs(:,[1 4 2 5 3 6]);
        end

        function k = computePairHessian(obj,X)
            h = 1e-7*max([1;abs(X(:))]);
            k = zeros(size(X,1),6,6);
            for j = 1:6
                e = zeros(1,6);
                e(j) = h;
                k(:,:,j) = (obj.computePairGradient(X+e) - obj.computePairGradient(X-e))/(2*h);
            end
            k = (k + permute(k,[1 3 2]))/2;
        end

        function [node,edge] = computeCandidatePairs(obj,x,margin)
            p = x(obj.edges(:,1),:);
            q = x(obj.edges(:,2),:);
            h = max(sqrt(sum((q-p).^2,2))) + margin;
            [edge,k] = obj.findNeighbours((p+q)/2,x(obj.boundaryNodes,:),h);
            node = obj.boundaryNodes(k);
        end

        function [i,j] = findNeighbours(obj,queries,points,h)
            pointKey = obj.computeCellKey(floor(points/h));
            [sortedKey,order] = sort(pointKey);
            [cellKeys,first]  = unique(sortedKey,'first');
            [~,last]          = unique(sortedKey,'last');
            [ox,oy] = meshgrid(-1:1,-1:1);
            nQueries = size(queries,1);
            queryOf  = repmat((1:nQueries)',9,1);
            queryKey = obj.computeCellKey(repmat(floor(queries/h),9,1) + repelem([ox(:),oy(:)],nQueries,1));
            [found,loc] = ismember(queryKey,cellKeys);
            queryOf = queryOf(found);
            loc     = loc(found);
            count   = last(loc) - first(loc) + 1;
            i       = repelem(queryOf,count);
            start   = repelem(first(loc),count);
            offset  = (1:sum(count))' - repelem(cumsum(count)-count,count) - 1;
            j       = order(start + offset);
        end

        function n = nDofs(obj)
            n = 2*obj.mesh.nnodes;
        end

    end

    methods (Access = private)

        function computeOrientedBoundaryEdges(obj)
            obj.mesh.computeEdges();
            e = obj.mesh.edges;
            nElemsOfEdge = accumarray(e.edgesInElem(:),1);
            [elem,local] = find(ismember(e.edgesInElem,find(nElemsOfEdge == 1)));
            ids   = e.edgesInElem(sub2ind(size(e.edgesInElem),elem,local));
            ab    = e.nodesInEdges(ids,:);
            other = sum(obj.mesh.connec(elem,:),2) - sum(ab,2);
            x     = obj.mesh.coord;
            isRight = obj.cross2D(x(ab(:,2),:)-x(ab(:,1),:),x(other,:)-x(ab(:,1),:)) < 0;
            ab(isRight,:) = ab(isRight,[2 1]);
            obj.edges = ab;
            obj.boundaryNodes = unique(ab(:));
        end

        function itDoes = doSegmentsCross(obj,p1,q1,p2,q2)
            d1 = q1 - p1;
            d2 = q2 - p2;
            r  = p2 - p1;
            den = obj.cross2D(d1,d2);
            a = obj.cross2D(r,d2)./den;
            b = obj.cross2D(r,d1)./den;
            itDoes = a > 0 & a < 1 & b > 0 & b < 1;
        end

    end

    methods (Static, Access = protected)

        function X = gatherPairCoordinates(x,node,a,b)
            X = [x(node,:),x(a,:),x(b,:)];
        end

        function c = cross2D(u,v)
            c = u(:,1).*v(:,2) - u(:,2).*v(:,1);
        end

        function key = computeCellKey(cell)
            key = (cell(:,1)+1e6)*2e6 + (cell(:,2)+1e6);
        end

    end

end
