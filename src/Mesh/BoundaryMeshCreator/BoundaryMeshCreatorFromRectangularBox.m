classdef BoundaryMeshCreatorFromRectangularBox < BoundaryMeshCreator

    properties (Access = private)
        nSides
        nFaces
        nDim
    end

    properties (Access = private)
        backgroundMesh
        dimension
        forUnfitted
    end

    methods (Access = public)

        function obj = BoundaryMeshCreatorFromRectangularBox(cParams)
            obj.init(cParams)
        end

        function b = create(obj)
            bMeshes = cell(obj.nFaces,1);
            for iDime = 1:obj.nDim
                for iSide = 1:obj.nSides
                    iFace = obj.computeIface(iSide,iDime);
                    bMeshes{iFace} = obj.createBoundaryMesh(iDime,iSide);
                end
            end
            b = bMeshes;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.dimension     = cParams.dimension;
            obj.nSides = 2;
            obj.nDim   = obj.backgroundMesh.ndim + obj.backgroundMesh.kFace;
            obj.nFaces = obj.nDim*obj.nSides;
            % obj.forUnfitted = cParams.forUnfitted;
        end

        function m = createBoundaryMesh(obj,iDime,iSide)
            nodes       = obj.obtainBoxNodes(iDime,iSide);
            coords      = obj.computeCoords(nodes);
            connec      = obj.computeConnectivities(nodes,iDime,obj.backgroundMesh.type);
            s.connec      = connec;
            s.coord       = coords;
            s.nodesInBoxFaces = nodes;
            s.dimension  = iDime;
            s.kFace      = obj.backgroundMesh.kFace;
            s.isRectangularBox = true;
            m = BoundaryMesh(s);
        end

        function coords = computeCoords(obj,nodes)
            coords = obj.backgroundMesh.coord(nodes,:);
        end

        function connec = computeConnectivities(obj,nodes,iDime,type)
            facetCoords = obj.computeFacetCoords(nodes,iDime);
            switch obj.nDim
                case 2
                    connec = obj.computeConnectivities1D(facetCoords);
                case 3
                    tf = all(facetCoords(:,1) == facetCoords(1,1)) || all(facetCoords(:,2) == facetCoords(1,2));
                    if tf
                        connec = obj.computeConnectivities1D(facetCoords);
                    else
                        connec = obj.computeConnectivities2D(facetCoords,type);
                    end
            end


        end

        function connec = computeConnecQuads(obj,coord)
            x = coord(:,1);
            y = coord(:,2);

            ux = unique(x);
            uy = unique(y);

            nx = numel(ux);
            ny = numel(uy);

            [~, ix] = ismember(x, ux);
            [~, iy] = ismember(y, uy);

            % Build index grid
            gridIndex = zeros(ny, nx);
            ind = sub2ind([ny nx], iy, ix);
            gridIndex(ind) = 1:length(x);

            A = gridIndex(1:end-1, 1:end-1);
            B = gridIndex(1:end-1, 2:end);
            C = gridIndex(2:end,   2:end);
            D = gridIndex(2:end,   1:end-1);

            connec = [A(:), B(:), C(:), D(:)];
            % % % Base indices (top-left corner of each quad)
            % % i = 1:nx-1;
            % % j = 1:ny-1;
            % %
            % % [I,J] = meshgrid(i,j);
            % % n1 = (J-1)*nx + I;
            % %
            % % % Build all quads at once
            % % quads = [ ...
            % %     n1(:), ...
            % %     n1(:)+1, ...
            % %     n1(:)+nx+1, ...
            % %     n1(:)+nx ];
        end

        function nodes = obtainBoxNodes(obj,iDime,iSide)
            dim = obj.dimension(iDime);
            coordDim = obj.backgroundMesh.coord(:,dim);
            switch iSide
                case 1
                    xL = min(coordDim);
                case 2
                    xL = max(coordDim);
            end
            % nodes = coordDim == xL;
            nodes = abs(coordDim - xL) <= 1e-11;
        end

        function facetCoord = computeFacetCoords(obj,nodes,idime)
            coord      = obj.backgroundMesh.coord(nodes,:);
            facetDim   = setdiff(1:obj.nDim,idime);
            facetCoord = coord(:,obj.dimension(facetDim));
        end

        function iFace = computeIface(obj,iSide,iDime)
            iFace = (iDime-1)*obj.nSides + iSide;
        end

        function connec = computeConnectivities2D(obj,coord,type)
            if strcmp('TETRAHEDRA',type)
                DT = delaunayTriangulation(coord);
                connec = DT.ConnectivityList;
            elseif strcmp('HEXAHEDRA',type)
                connec = obj.computeConnecQuads(coord);
            end
        end

    end

    methods (Access = private, Static)

        function connec = computeConnectivities1D(coord)
            [~,I] = sort(coord);
            connec = [I circshift(I,-1)];
            connec(end,:) = [];
        end

    end

end