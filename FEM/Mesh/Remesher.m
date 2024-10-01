classdef Remesher < handle

    properties (Access = public)
        %    coord
        %    connec
        fineMesh
        fineContMesh
    end

    properties (Access = private)
        mesh
        cellsToRemesh
    end

    methods (Access = public)

        function obj = Remesher(cParams)
            obj.init(cParams)
        end

        function m = remesh(obj)
            s.coord  = obj.computeCoords();
            s.connec = obj.computeConnectivities();
            m = Mesh.create(s);            
        end

    end

    methods (Access = public, Static)

        function mF = compute(m)
            s.mesh = m;
            r = Remesher(s);
            mF = r.remesh();
        end

    end
   
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.cellsToRemesh = 1:obj.mesh.nelem;
        end

        function coord = computeCoords(obj)
            oldCoord = obj.mesh.coord;
            newCoord = obj.computeNewCoord(obj.mesh,obj.mesh.coord);
            coord = [oldCoord;newCoord];
        end

        function f = computeNewCoord(obj,m,fNodes)
            s.edgeMesh = m.computeEdgeMesh();
            s.fNodes   = fNodes;
            eF         = EdgeFunctionInterpolator(s);
            f = eF.compute()';
        end

        function allConnec = computeConnectivities(obj)
            oldConnec = obj.mesh.connec;
            newConnec = obj.computeNewConnectivities();
            allConnec = [oldConnec;newConnec];
        end

       function newConnec = computeNewConnectivities(obj)
            [nV1,nV2,nV3] = obj.computeNodeInVertex();
            [nE1,nE2,nE3] = obj.computeNodeInEdge();
            connec(:,1,:) = [nV1 nE1 nE3];
            connec(:,2,:) = [nE1 nV2 nE2];
            connec(:,3,:) = [nE1 nE2 nE3];
            connec(:,4,:) = [nE3 nE2 nV3];
            newConnec = reshape(connec,[],3);
        end

        function [nV1,nV2,nV3] = computeNodeInVertex(obj)
            cells        = obj.cellsToRemesh;
            vertexInCell = obj.mesh.connec(cells,:);
            nV1 = vertexInCell(:,1);
            nV2 = vertexInCell(:,2);
            nV3 = vertexInCell(:,3);
        end

        function [nE1,nE2,nE3] = computeNodeInEdge(obj)
            e  = obj.computeEdgesInCell();
            e1 = e(:,1);
            e2 = e(:,2);
            e3 = e(:,3);
            newNodes = obj.computeNewNodes();
            nE1 = newNodes(e1);
            nE2 = newNodes(e2);
            nE3 = newNodes(e3);
        end

        function eC = computeEdgesInCell(obj)
            obj.mesh.computeEdges();
            e = obj.mesh.edges;
            cells = obj.cellsToRemesh;
            eC =  e.edgesInElem(cells,:);
        end

        function newNodes = computeNewNodes(obj)
            e      = obj.mesh.edges;
            nEdges = size(e.nodesInEdges,1);
            nnodes        = obj.mesh.nnodes;
            newNodes(:,1) = (nnodes+1):(nnodes+nEdges);
        end

    end

end