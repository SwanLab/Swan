classdef Mesh < handle

    properties (GetAccess = public, SetAccess = private)
        type
        kFace

        coord
        connec

        ndim
        nelem
        nnodes
        nnodeElem

       % remove (xFE)
        interpolation

        edges
        faces
        boundaryNodes
        boundaryElements
        coordElem
    end

    properties (Access = private)
        xVOld
        dVOld
        
    end    

    properties (Access = protected)
        xFE
    end

    methods (Static, Access = public)
        
        function obj = create(cParams)
            if ~(isfield(cParams,'kFace'))
                cParams.kFace = 0;
            end
            g = GeometryType.compute(cParams);
            cParams.type = MeshTypeComputer.compute(cParams.connec,g);
            switch g
                case 'Line'
                    obj = LineMesh(cParams);
                case 'Surface'
                    obj = SurfaceMesh(cParams);
                case 'Volume'
                    obj = VolumeMesh(cParams);
            end
        end
        
        function obj = createFromGiD(filename)
            reader = FemInputReaderGiD();
            a = reader.read(filename);
            obj = a.mesh;
        end

    end

    methods (Access = public)

        function L = computeCharacteristicLength(obj)
            xmin = min(obj.coord);
            xmax = max(obj.coord);
            L = norm(xmax-xmin);
        end

        function xV = computeBaricenter(obj)
            q = Quadrature.create(obj, 0);
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
            quad = Quadrature.create(obj, 1);
            volume = obj.computeDvolume(quad)';
            L = obj.computeSquarePerimeter();
            q = 4*sqrt(3)*volume./L;
        end

        function v = computeVolume(obj) % computeMeasure
            quad = Quadrature.create(obj,2);
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
                    n(:,1:obj.nelem) = squeeze(normal);
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

        function plotAllNodes(obj)
            nodes = 1:obj.nnodes;
            obj.plotNodes(nodes,'blue')
        end

        function plotNodes(obj,indeces,colorValue)
            ind(:,1) = indeces;
            b = num2str(ind);
            c = cellstr(b);
            dx = 0.01; dy = 0.01;
            x = obj.coord(ind,1)';
            y = obj.coord(ind,2)';
            t = text(x+dx,y+dy,c);
            set(t,'Color',colorValue)
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

        function mF = remesh(obj) % only tri mesh
            % for quad, QuadToTriMeshConverter
            mC = obj;
            mF = Remesher.compute(mC);
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

        function dV = computeDvolume(obj,quad)
            xV = quad.posgp;
            if ~isequal(xV,obj.xVOld)
                w = reshape(quad.weigp,[quad.ngaus 1]);
                dVolume = w.*obj.computeJacobianDeterminant(quad.posgp);
                dV = reshape(dVolume, [quad.ngaus, obj.nelem]);
                obj.dVOld = dV;
                obj.xVOld = xV;
            else
                dV = obj.dVOld;
            end
        end

        %% Remove

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
                    validNodes = find(domain(obj.coord));
                    validElems = find(sum(ismember(obj.connec, validNodes),2) == obj.nnodeElem); % == 2 because line
                    coord_valid  = obj.coord(validNodes, :);
                    connec_valid = obj.connec(validElems,:);
                    % connecGlobal = l2gBound(connec_valid);
                    
                    newNodes = (1:size(coord_valid,1))';
                    boundary2local(validNodes) = newNodes;
                    newConnec = boundary2local(connec_valid);
                    
                    s.connec = newConnec;
                    s.coord = coord_valid;
                    s.kFace = -1;
                    
                    m = Mesh.create(s);
                    l2g(newNodes(:)) = validNodes;

%                     error('Cannot yet get boundary submesh for 3D')
            end
        end

    end

    methods (Access = protected)
        function obj = Mesh(cParams)
            obj.init(cParams);
            obj.computeDimensionParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();
        end
    end

    methods (Access = public) % ?????????

        function J = computeJacobian(obj,xV)
            nDimGlo  = size(obj.coordElem,1);
            nElem    = size(obj.coordElem,3);
            dShapes  = obj.interpolation.computeShapeDerivatives(xV);
            nDimElem = size(dShapes,1);
            nPoints  = size(xV,2);
            J = zeros(nDimElem,nDimGlo,nPoints,nElem);
            for iDimGlo = 1:nDimGlo
                for iDimElem = 1:nDimElem
                    dShapeIK = squeezeParticular(dShapes(iDimElem,:,:),1)';
                    xKJ = squeezeParticular(obj.coordElem(iDimGlo,:,:),1);
                    jacIJ    = dShapeIK*xKJ;
                    J(iDimElem,iDimGlo,:,:) = squeezeParticular(J(iDimElem,iDimGlo,:,:),[1 2]) + jacIJ;
                end
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

        function createInterpolation(obj)
            obj.interpolation = Interpolation.create(obj.type,'LINEAR');
        end

        function computeElementCoordinates(obj)
            obj.computeCoordFEfunction();
            obj.coordElem = obj.xFE.getFvaluesByElem();
        end

        function computeCoordFEfunction(obj)
            s.mesh    = obj;
            s.order   = 'P1';
            s.fValues = obj.coord;
            coordP1   = LagrangianFunction(s);
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
