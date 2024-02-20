classdef Mesh < handle

    properties (GetAccess = public, SetAccess = private)
        type
        kFace
        % geometryType

        coord
        connec

        ndim
        nelem
        nnodes
        nnodeElem

        coordElem % remove (xFE)
        interpolation

        edges
        faces
        boundaryNodes
        boundaryElements
    end

    properties (Access = private)
        xFE
        geometry
    end

    methods (Static, Access = public)
        
        function obj = create(cParams)
            s = SettingsMesh(cParams);
            switch s.geometryType
                case 'Line'
                    obj = LineMesh(s);
                case 'Surface'
                    obj = SurfaceMesh(s);
                case 'Volume'
                    obj = VolumeMesh(s);
            end
        end

    end

    methods (Access = public)

        function obj = Mesh(cParams)
            obj.init(cParams);
            obj.computeDimensionParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();
            obj.createGeometry();
        end

        function L = computeCharacteristicLength(obj)
            xmin = min(obj.coord);
            xmax = max(obj.coord);
            L = norm(xmax-xmin);
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
        
        function hMin = computeMinCellSize(obj)
            if isempty(obj.edges)
                obj.computeEdges();
            end
            x1 = obj.coord(obj.edges.nodesInEdges(:,1),:);
            x2 = obj.coord(obj.edges.nodesInEdges(:,2),:);
            x1x2 = (x2-x1);
            hMin = min(sqrt(sum(x1x2.^2,2)));
        end

        function hMean = computeMeanCellSize(obj)
            if isempty(obj.edges)
                obj.computeEdges();
            end
            x1 = obj.coord(obj.edges.nodesInEdges(:,1),:);
            x2 = obj.coord(obj.edges.nodesInEdges(:,2),:);
            x1x2 = (x2-x1);
            hMean = mean(sqrt(sum(x1x2.^2,2)));
        end

        function hMax = computeMaxCellSize(obj)
            if isempty(obj.edges)
                obj.computeEdges();
            end
            x1 = obj.coord(obj.edges.nodesInEdges(:,1),:);
            x2 = obj.coord(obj.edges.nodesInEdges(:,2),:);
            x1x2 = (x2-x1);
            hMax = max(sqrt(sum(x1x2.^2,2)));
        end

        function q = computeElementQuality(obj) % check for 3d
            quad = Quadrature.set(obj.type);
            quad.computeQuadrature('CONSTANT');
            volume = obj.computeDvolume(quad);
            L(1,:) = obj.computeSquarePerimeter();
            q = 4*sqrt(3)*volume./L;
        end

        function v = computeVolume(obj) % computeMeasure
            quad = Quadrature.set(obj.type);
            quad.computeQuadrature('CONSTANT');
            v = obj.computeDvolume(quad);
            v = sum(v(:));
        end

        function computeEdges(obj) % nonsense for lines
            s.nodesByElem = obj.connec;
            s.type = obj.type;
            edge = EdgesConnectivitiesComputer(s);
            edge.compute();
            obj.edges = edge;
        end
        
        function computeFaces(obj) % nonsense for lines
            s.nodesByElem = obj.connec;
            s.type = obj.type;
            face = FacesConnectivitiesComputer(s);
            face.compute();
            obj.faces = face;
        end

        function eM = computeEdgeMesh(obj) % nonsense for lines
            obj.computeEdges();
            s.coord  = obj.coord;
            s.connec = obj.edges.nodesInEdges;
            s.kFace  = obj.kFace -1;
            eM = Mesh.create(s);
        end

        function m = computeCanonicalMesh(obj)
            s.remainingNodes = unique(obj.connec);
            s.mesh        = obj;
            c = CannonicalMeshComputer(s);
            m = c.compute();
        end

        function plotNormals(obj) % volume
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

        function cells = computeAllCellsOfVertex(obj,vertex)
            vertexInCell  = obj.connec;
            isInCell      = any(vertexInCell == vertex,2);
            allCells(:,1) = 1:size(isInCell,1);
            cells         = allCells(isInCell);
        end

        function cV = computeConnectedVertex(obj,vertex)
            cV  = obj.edges.computeConnectedVertex(vertex);
        end

        function m = remesh(obj,nLevels) % only tri mesh
            % for quad, QuadToTriMeshConverter
            s.mesh = obj;
            s.nLevels = nLevels;
            r = Remesher(s);
            m = r.compute();
        end

        function exportSTL(obj)
            s.mesh = obj;
            me = STLExporter(s);
            me.export();
        end

        function print(obj, filename, software)
            if nargin == 2; software = 'Paraview'; end
            p1 = LagrangianFunction.create(obj,1, 'P1');
            p1.print(filename, software);
        end

        %% Heavy refactoring

        % Separate Mesh into LineMesh, SurfaceMesh, VolumeMesh
        % DELETE Geometry

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

        function n = getNormals(obj) % only 
            quad = Quadrature.set(obj.type);
            quad.computeQuadrature('CONSTANT');
            g = obj.geometry;
            g.computeGeometry(quad,obj.interpolation);
            n = g.normalVector;
        end

        %% Remove

        function setCoord(obj,newCoord)
            obj.coord = newCoord;
        end

        function mD = createDiscontinuousMesh(obj) % P1D
            ndims = size(obj.coord, 2);
            nNodesDisc = obj.nnodeElem*obj.nelem;
            nodesDisc  = 1:nNodesDisc;
            connecDisc = reshape(nodesDisc,obj.nnodeElem,obj.nelem)';
            coordD = reshape(obj.xFE.fValues, [ndims, nNodesDisc])';
            s.connec = connecDisc;
            s.coord  = coordD;
            mD = Mesh.create(s);
        end

        function [m, l2g] = createSingleBoundaryMesh(obj)
            % To BoundaryMesh
            x = obj.coord(:,1);
            y = obj.coord(:,2);
            
            k = boundary(x,y);
            k = k(1:end-1);
            originalNodes = k;
            newNodes = (1:length(k))';
            boundaryCoords = [x(k), y(k)];
            boundaryConnec = [newNodes, circshift(newNodes,-1)];

            s.connec = boundaryConnec;
            s.coord = boundaryCoords;
            s.kFace = -1;
            
            m = Mesh.create(s);
            l2g(newNodes(:)) = originalNodes(:);
        end
        
        function [m, l2g] = getBoundarySubmesh(obj, domain)
            % To BoundaryMesh -- obj.boundary.mesh{iMesh} instead of
            % obj.boundary{iMesh}.mesh
            switch obj.ndim
                case 2
                    [mBound, l2gBound] = obj.createSingleBoundaryMesh();
                    validNodes = find(domain(mBound.coord));
                    validElems = find(sum(ismember(mBound.connec, validNodes),2) == 2); % == 2 because line
                    coord_valid  = mBound.coord(validNodes, :);
                    connec_valid = mBound.connec(validElems,:);
                    % connecGlobal = l2gBound(connec_valid);
                    
                    newNodes = (1:size(coord_valid,1))';
                    boundary2local(validNodes) = newNodes;
                    newConnec = boundary2local(connec_valid);
                    
                    s.connec = newConnec;
                    s.coord = coord_valid;
                    s.kFace = -1;
                    
                    m = Mesh.create(s);
                    l2g(newNodes(:)) = l2gBound(validNodes);
                otherwise
                    error('Cannot yet get boundary submesh for 3D')
            end
        end

    end

    methods (Access = private)

        function init(obj,s)
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

        function createGeometry(obj)
            s.mesh = obj;
            obj.geometry = Geometry.create(s);
        end

        function createInterpolation(obj)
            obj.interpolation = Interpolation.create(obj.type,'LINEAR');
        end

        function computeElementCoordinates(obj)
            obj.computeCoordFEfunction();
            obj.coordElem = obj.xFE.fValues;
        end

        function computeCoordFEfunction(obj)
            s.mesh    = obj;
            s.order   = 'P1';
            s.fValues = obj.coord;
            coordP1 = LagrangianFunction(s);
            obj.xFE = coordP1.project('P1D');
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

    end

end
