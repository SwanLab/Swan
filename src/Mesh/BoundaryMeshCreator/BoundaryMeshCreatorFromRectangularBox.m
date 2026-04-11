classdef BoundaryMeshCreatorFromRectangularBox < BoundaryMeshCreator
    
    properties (Access = private)
        nSides
        nFaces
        nDim
    end
    
    properties (Access = private)
        backgroundMesh
        dimension
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
        end
        
        function m = createBoundaryMesh(obj,iDime,iSide) 
            nodes       = obj.obtainBoxNodes(iDime,iSide);
            coords      = obj.computeCoords(nodes);
            connec      = obj.computeConnectivities(nodes,iDime);
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
        
        function connec = computeConnectivities(obj,nodes,iDime)
            facetCoords = obj.computeFacetCoords(nodes,iDime);
            switch obj.nDim
                case 2
                    connec = obj.computeConnectivities1D(facetCoords);
                case 3
                    % connec = obj.computeConnectivities2D(facetCoords);
                    tf = all(facetCoords(:,1) == facetCoords(1,1)) || all(facetCoords(:,2) == facetCoords(1,2));
                    if tf
                        connec = obj.computeConnectivities1D(facetCoords);
                    else
                        connec = obj.computeConnectivities2D(facetCoords);
                    end
            end
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
            nodes = coordDim == xL;
        end
        
        function facetCoord = computeFacetCoords(obj,nodes,idime)
            coord      = obj.backgroundMesh.coord(nodes,:);
            facetDim   = setdiff(1:obj.nDim,idime);
            facetCoord = coord(:,obj.dimension(facetDim));
        end

        function iFace = computeIface(obj,iSide,iDime)
            iFace = (iDime-1)*obj.nSides + iSide;
        end
        
    end
    
    methods (Access = private, Static)
        
        function connec = computeConnectivities1D(coord)
            [~,I] = sort(coord);
            connec = [I circshift(I,-1)];
            connec(end,:) = [];
        end
        
        function connec = computeConnectivities2D(coord)
            DT = delaunayTriangulation(coord);
            connec = DT.ConnectivityList;
        end
        
    end
    
end