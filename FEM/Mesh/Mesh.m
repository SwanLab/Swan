classdef Mesh < handle

    properties (GetAccess = public, SetAccess = private)
        type
        kFace
        geometryType

        coord
        connec

        ndim
        nelem
        nnodes
        nnodeElem

        coordElem
        interpolation

        edges
        boundaryNodes
        boundaryElements

        masterSlaveNodes
    end

    properties (Access = private)
        xFE
        geometry
        triMesh
    end


    methods (Access = public)

        function obj = Mesh(cParams)
            obj.init(cParams);
            obj.computeDimensionParams();
            obj.computeGeometryType();
            obj.computeType();
            obj.createInterpolation();
            obj.computeElementCoordinates();
            obj.createGeometry();
        end

        function obj = createFromFile(obj,cParams)
            testName = cParams.testName;
            [coordV, connecV] = obj.readCoordConnec(testName);
            s.coord  = coordV(:,2:end-1);
            s.connec = connecV(:,2:end);
            obj = Mesh(s);
        end

        function plot(obj)
            s.mesh = obj;
            s = SettingsMeshPlotter(s);
            mP = MeshPlotter(s);
            mP.plot();
        end

        function hMin = computeMinCellSize(obj)
            x1(:,1) = obj.coord(obj.connec(:,1),1);
            x1(:,2) = obj.coord(obj.connec(:,1),2);
            x2(:,1) = obj.coord(obj.connec(:,2),1);
            x2(:,2) = obj.coord(obj.connec(:,2),2);
            x3(:,1) = obj.coord(obj.connec(:,3),1);
            x3(:,2) = obj.coord(obj.connec(:,3),2);
            x1x2 = (x2-x1);
            x2x3 = (x3-x2);
            x1x3 = (x1-x3);
            n12 = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
            n23 = sqrt(x2x3(:,1).^2 + x2x3(:,2).^2);
            n13 = sqrt(x1x3(:,1).^2 + x1x3(:,2).^2);
            hs = min([n12,n23,n13],[],2);
            hMin = min(hs);
        end

        function hMean = computeMeanCellSize(obj)
            switch obj.type
                case {'LINE'}
                    x1(:,1) = obj.coord(obj.connec(:,1),1);
                    x1(:,2) = obj.coord(obj.connec(:,1),2);
                    x2(:,1) = obj.coord(obj.connec(:,2),1);
                    x2(:,2) = obj.coord(obj.connec(:,2),2);
                    x1x2 = (x2-x1);                    
                    hs = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
                    hMean = max(hs);  
                otherwise
                    x1(:,1) = obj.coord(obj.connec(:,1),1);
                    x1(:,2) = obj.coord(obj.connec(:,1),2);
                    x2(:,1) = obj.coord(obj.connec(:,2),1);
                    x2(:,2) = obj.coord(obj.connec(:,2),2);
                    x3(:,1) = obj.coord(obj.connec(:,3),1);
                    x3(:,2) = obj.coord(obj.connec(:,3),2);
                    x1x2 = (x2-x1);
                    x2x3 = (x3-x2);
                    x1x3 = (x1-x3);
                    n12 = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
                    n23 = sqrt(x2x3(:,1).^2 + x2x3(:,2).^2);
                    n13 = sqrt(x1x3(:,1).^2 + x1x3(:,2).^2);
                    hs = max([n12,n23,n13],[],2);
                    hMean = max(hs);                  
            end
        end

        function L = computeCharacteristicLength(obj)
            xmin = min(obj.coord);
            xmax = max(obj.coord);
            L = norm(xmax-xmin);
        end

        function setCoord(obj,newCoord)
            obj.coord = newCoord;
        end

        function setConnec(obj,newConnec)
            obj.connec = newConnec;
        end

        function xV = computeBaricenter(obj)
            q = Quadrature.set(obj.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            xV = squeeze(obj.xFE.evaluate(xV));
        end

        function xGauss = computeXgauss(obj,xV)
            xGauss = obj.xFE.evaluate(xV);
        end

        function dvolume = computeDvolume(obj,quad)
            g = obj.geometry;
            g.computeGeometry(quad,obj.interpolation);
            dvolume = g.dvolu;
            dvolume = dvolume';
        end

        function invJac = computeInverseJacobian(obj,quad,int)
            g = obj.geometry;
            invJac = g.computeInverseJacobian(quad,int);
        end

        function n = getNormals(obj)
            quad = Quadrature.set(obj.type);
            quad.computeQuadrature('CONSTANT');
            g = obj.geometry;
            g.computeGeometry(quad,obj.interpolation);
            n = g.normalVector;
        end

        function q = computeElementQuality(obj)
            quad = Quadrature.set(obj.type);
            quad.computeQuadrature('CONSTANT');
            volume = obj.computeDvolume(quad);
            L(1,:) = obj.computeSquarePerimeter();
            q = 4*sqrt(3)*volume./L;
        end

        function v = computeVolume(obj)
            quad = Quadrature.set(obj.type);
            quad.computeQuadrature('CONSTANT');
            v = obj.computeDvolume(quad);
            v = sum(v(:));
        end

        function computeEdges(obj)
            s.nodesByElem = obj.connec;
            edge = EdgesConnectivitiesComputer(s);
            edge.compute();
            obj.edges = edge;
        end

        function eM = computeEdgeMesh(obj)
            obj.computeEdges;
            s.coord  = obj.coord;
            s.connec = obj.edges.nodesInEdges;
            s.kFace  = obj.kFace -1;
            eM = Mesh(s);
        end

        function m = computeCanonicalMesh(obj)
            s.remainingNodes = unique(obj.connec);
            s.mesh        = obj;
            c = CannonicalMeshComputer(s);
            m = c.compute();
        end

        function setMasterSlaveNodes(obj,nodes)
            obj.masterSlaveNodes = nodes;
        end

        function computeMasterSlaveNodes(obj)
           mR = MasterSlaveRelator(obj.coord);
           nodes = mR.getRelation();
           obj.masterSlaveNodes = nodes;
        end

        function plotNormals(obj)
            switch obj.ndim
                case 3
                    normal = obj.getNormals();
                    n(:,1:obj.nelem) = squeeze(normal)';
                    xy = obj.computeBaricenter();
                    q = quiver3(xy(1,:),xy(2,:),xy(3,:),n(1,:),n(2,:),n(3,:),'k');
                    h = obj.computeMeanCellSize();
                    q.LineWidth = h;
                    q.Marker = '.';
                    q.MaxHeadSize = 1;
                    ah = annotation('arrow','headStyle','cback1');
                    set(ah,'parent',gca);
                    axis equal
            end
        end

        function bMesh = createBoundaryMesh(obj)
            if isempty(obj.boundaryNodes) || isempty(obj.boundaryElements)
                s.backgroundMesh = obj;
                s.dimension = 1:obj.ndim;
                s.type = 'FromReactangularBox';
                bC = BoundaryMeshCreator.create(s);
                bMesh = bC.create();
            else
                s.borderNodes    = obj.boundaryNodes;
                s.borderElements = obj.boundaryElements;
                s.backgroundMesh = obj;
                s.type = 'FromData';
                b = BoundaryMeshCreator.create(s);
                bMesh = b.create();
            end
        end

        function mD = createDiscontinuousMesh(obj)
            ndims = size(obj.coord, 2);
            nNodesDisc = obj.nnodeElem*obj.nelem;
            nodesDisc  = 1:nNodesDisc;
            connecDisc = reshape(nodesDisc,obj.nnodeElem,obj.nelem)';
            coordD = reshape(obj.xFE.fValues, [ndims, nNodesDisc])';
            s.connec = connecDisc;
            s.coord  = coordD;
            mD = Mesh(s);
        end

        function cells = computeAllCellsOfVertex(obj,vertex)
            vertexInCell  = obj.connec;
            isInCell      = any(vertexInCell == vertex,2);
            allCells(:,1) = 1:size(isInCell,1);
            cells         = allCells(isInCell);
        end

        function cV = computeConnectedVertex(obj,vertex)
            cV  = obj.edges.computeConnectedVertex(vertex);
        end

        function m = remesh(obj,nLevels)
            s.mesh = obj;
            s.nLevels = nLevels;
            r = Remesher(s);
            m = r.compute();
        end

        function exportSTL(obj, file)
            obj.triangulateMesh();
            stlwrite(obj.triMesh, [file '.stl'])
        end

        function exportSTL3D(obj, file)
            model = createpde;
            g     = importGeometry(model,[file,'.stl']);
            g_3D  = extrude(g,5);
            pdegplot(g_3D,"FaceLabels","on")
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            s = SettingsMesh(cParams);
            if isfield(cParams,'boundaryNodes')
               obj.boundaryNodes = cParams.boundaryNodes;
            end
            if isfield(cParams,'boundaryElements')
               obj.boundaryElements = cParams.boundaryElements;
            end
            if isfield(cParams,'masterSlaveNodes')
               obj.masterSlaveNodes = cParams.masterSlaveNodes;
            end
            obj.coord  = s.coord;
            obj.connec = s.connec;
            obj.type   = s.type;
            obj.kFace  = s.kFace;
        end

        function computeDimensionParams(obj)
            obj.nnodes = size(obj.coord,1);
            obj.ndim  = size(obj.coord,2);
            obj.nelem = size(obj.connec,1);
            obj.nnodeElem = size(obj.connec,2);
        end

        function computeType(obj)
            s.geometryType = obj.geometryType;
            s.nnodeElem    = obj.nnodeElem;
            t = MeshTypeComputer(s);
            obj.type = t.compute();
        end

        function computeGeometryType(obj)
            s.ndim  = obj.ndim;
            s.kFace = obj.kFace;
            g = GeometryTypeComputer(s);
            obj.geometryType = g.compute();
        end

        function createGeometry(obj)
            s.mesh = obj;
            obj.geometry = Geometry.create(s);
        end

        function createInterpolation(obj)
            obj.interpolation = Interpolation.create(obj,'LINEAR');
        end

        function computeElementCoordinates(obj)
            obj.computeCoordFEfunction();
            obj.coordElem = obj.xFE.fValues;
        end

        function computeCoordFEfunction(obj)
            s.mesh    = obj;
            s.fValues = obj.coord;
            coordP1 = P1Function(s);
            obj.xFE = obj.projectToP1Discontinuous(coordP1);
        end

        function p1d = projectToP1Discontinuous(obj, f)
            s.mesh   = obj;
            sP.origin = 'P1';
            sP.x = f;
            p = ProjectorToP1discont(s);
            p1d = p.project(sP);
        end

        function L = computeSquarePerimeter(obj)
            obj.computeEdges();
            nElem = size(obj.connec,1);
            L = zeros(nElem,1);
            for iedge = 1:obj.edges.nEdgeByElem
                edge = obj.edges.edgesInElem(:,iedge);
                nodesEdge = obj.edges.nodesInEdges(edge,:);
                for idim = 1:obj.ndim
                    xA = obj.coord(nodesEdge(:,1),idim);
                    xB = obj.coord(nodesEdge(:,2),idim);
                    L = L + (xA - xB).^2;
                end
            end
        end

        function triangulateMesh(obj)
            P = obj.coord;
            T = obj.connec;
            obj.triMesh = triangulation(T,P);
        end

    end

    methods (Access = private, Static)

        function [coord, connec] = readCoordConnec(testName)
            run(testName)
        end

    end

end
