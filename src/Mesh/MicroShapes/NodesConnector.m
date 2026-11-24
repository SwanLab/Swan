classdef NodesConnector < handle

    properties (Access = private)
        nodes
        coord
        div
        T
    end

    properties (Constant, Access = private)
        orthogonalityTolerance = 1e-12
    end

    properties (Access = public)
        connec
    end

    methods (Access = public)

        function obj = NodesConnector(cParams)
            obj.init(cParams);
        end

        function computeConnections(obj)
            if obj.nodes.vert == 4
                obj.connec = obj.computeStructuredConnections();
            else
                obj.connec = delaunay(obj.coord);
                obj.deleteBoundaryOnlyElements();
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nodes = cParams.nodes;
            obj.coord = cParams.coord;
            obj.div   = cParams.div;
            obj.T     = eye(2);
            if isfield(cParams,'T')
                obj.T = cParams.T;
            end
        end

        function connec = computeStructuredConnections(obj)
            n       = obj.div(:)';
            nodeId  = obj.computeNodeGrid(n);
            [i,j]   = ndgrid(1:n(1),1:n(2));
            n00     = nodeId(sub2ind(n+1,i,j));
            n10     = nodeId(sub2ind(n+1,i+1,j));
            n11     = nodeId(sub2ind(n+1,i+1,j+1));
            n01     = nodeId(sub2ind(n+1,i,j+1));
            useMain = obj.selectMainDiagonal(i(:),j(:));
            mainA   = [n00(:),n10(:),n11(:)];
            mainB   = [n00(:),n11(:),n01(:)];
            antiA   = [n00(:),n10(:),n01(:)];
            antiB   = [n10(:),n11(:),n01(:)];
            connec  = [mainA(useMain,:);mainB(useMain,:);antiA(~useMain,:);antiB(~useMain,:)];
        end

        function useMain = selectMainDiagonal(obj,i,j)
            [e1,e2]  = obj.computePhysicalEdges();
            cosAngle = dot(e1,e2)/(norm(e1)*norm(e2));
            if abs(cosAngle) < obj.orthogonalityTolerance
                useMain = mod(i+j,2) == 0;
            else
                useMain = repmat(cosAngle < 0,size(i));
            end
        end

        function [e1,e2] = computePhysicalEdges(obj)
            origin = obj.coord(1,:);
            e1     = obj.T*(obj.coord(2,:)-origin)';
            e2     = obj.T*(obj.coord(4,:)-origin)';
        end

        function nodeId = computeNodeGrid(obj,n)
            gridIndex = obj.computeGridIndex(n);
            nodeId    = zeros(n+1);
            position  = sub2ind(n+1,gridIndex(:,1)+1,gridIndex(:,2)+1);
            nodeId(position) = 1:size(obj.coord,1);
            if any(nodeId(:) == 0) || numel(unique(position)) ~= size(obj.coord,1)
                error('NodesConnector:InvalidGrid','Nodes do not form a structured grid.');
            end
        end

        function gridIndex = computeGridIndex(obj,n)
            origin    = obj.coord(1,:);
            basis     = [obj.coord(2,:)-origin;obj.coord(4,:)-origin]';
            xi        = (basis\(obj.coord-origin)')';
            gridIndex = round(xi.*n);
        end

        function deleteBoundaryOnlyElements(obj)
            isBoundary = obj.connec <= obj.nodes.bound;
            obj.connec(all(isBoundary,2),:) = [];
        end

    end

end