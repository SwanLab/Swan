classdef BoundaryMeshCreatorFromRectangularBox < handle
    
    properties (Access = private)
        nSides
        nFaces
    end
    
    properties (Access = private)
        backgroundMesh
    end
    
    methods (Access = public)
        
        function obj = BoundaryMeshCreatorFromRectangularBox(cParams)
            obj.init(cParams)
        end
        
        function b = create(obj)
            bMeshes = cell(obj.nFaces,1);
            for iDime = 1:obj.backgroundMesh.ndim
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
            obj.nSides = 2;
            obj.nFaces = obj.backgroundMesh.ndim*obj.nSides;
        end
        
        function m = createBoundaryMesh(obj,iDime,iSide)            
            nodes  = obj.obtainBoxNodes(iDime,iSide);
            coords = obj.computeCoords(nodes);            
            connec = obj.computeConnectivities(coords,iDime);
            s.connec = connec;
            s.coord  = coords;
            s.nodesInBoxFaces = nodes;
            m = BoundaryMesh(s);
        end       
        
        function coords = computeCoords(obj,nodes)
            coords = obj.backgroundMesh.coord(nodes,:);            
        end
        
        function connec = computeConnectivities(obj,boxFaceCoords,iDime)
            facetCoords = obj.computeFacetCoords(boxFaceCoords,iDime);            
            switch obj.backgroundMesh.ndim
                case 2
                    connec = obj.computeConnectivities1D(facetCoords);
                case 3
                    connec = obj.computeConnectivities2D(facetCoords);
            end            
        end
        
        function nodes = obtainBoxNodes(obj,iDime,iSide)
            coordDim = obj.backgroundMesh.coord(:,iDime);
            switch iSide 
                case 1
                    xL = min(coordDim);
                case 2
                    xL = max(coordDim);
            end
            nodes = coordDim == xL;
        end
        
        function facetCoord = computeFacetCoords(obj,facetCoord,idime)
            nDim       = obj.backgroundMesh.ndim;
            facetDim   = setdiff(1:nDim,idime);
            facetCoord = facetCoord(:,facetDim);
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