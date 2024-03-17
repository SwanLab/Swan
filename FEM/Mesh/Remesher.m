classdef Remesher < handle

    properties (Access = public)
        %    coord
        %    connec
        fineMesh
    end

    properties (Access = private)
        mesh
        nLevels
        cellsToRemesh
    end

    methods (Access = public)

        function obj = Remesher(cParams)
            obj.init(cParams)
        end

        function mF = compute(obj)
            s.coord  = obj.computeCoords();
            s.connec = obj.computeConnectivities();
            mF = Mesh.create(s);
        end

        function remesh(obj)
            %  s.nodesByElem = obj.connec;
            %  edge = EdgesConnectivitiesComputer(s);
            %  edge.compute();
            m0 = obj.mesh;
            for iLevel = 1:obj.nLevels
                s.coord  = obj.computeCoords();
                s.connec = obj.computeConnectivities();
                m = Mesh.create(s);
                m = m.createDiscontinuousMesh();
                obj.mesh = m;
                obj.cellsToRemesh = 1:obj.mesh.nelem;
            end
            obj.fineMesh = obj.mesh;
            obj.mesh     = m0;
            obj.cellsToRemesh = 1:obj.mesh.nelem;
        end

        function f = interpolate(obj,f)
            m0 = obj.mesh;
            for iLevel = 1:obj.nLevels
                m = obj.mesh;
                s.coord  = obj.computeCoords();
                s.connec = obj.computeConnectivities();
                mF = Mesh.create(s);
                f  = f.refine(m,mF);
                mD = mF.createDiscontinuousMesh();
                obj.mesh = mD;
                obj.cellsToRemesh = 1:obj.mesh.nelem;
            end
            obj.mesh     = m0;
            obj.cellsToRemesh = 1:obj.mesh.nelem;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            %    obj.coord    = cParams.mesh.coord;
            %  obj.connec   = cParams.mesh.connec;
            obj.mesh     = cParams.mesh;
            obj.nLevels  = cParams.nLevels;
            obj.cellsToRemesh = 1:obj.mesh.nelem;
        end

        function coord = computeCoords(obj)
            oldCoord = obj.mesh.coord;
            newCoord = obj.computeNewCoord();
            coord = [oldCoord;newCoord];
        end

        function nC = computeNewCoord(obj)
            me           = obj.mesh.computeEdgeMesh();
            coordInEdges = me.computeBaricenter()';
            nC           = coordInEdges;
        end

        function allConnec = computeConnectivities(obj)
            oldConnec = obj.computeOldConnectivities();
            newConnec = obj.computeNewConnectivities();
            allConnec = [oldConnec,newConnec];
        end

        function oldConnec = computeOldConnectivities(obj)
            cells         = obj.cellsToRemesh;
            cellToNotMesh = setdiff(1:obj.mesh.nelem,cells);
            oldConnec     = obj.mesh.connec(cellToNotMesh,:);
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