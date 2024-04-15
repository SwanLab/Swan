classdef StructuredMesh < handle
    
    properties (Access = public)
        nDivs
        divs
        mesh
    end
    
    properties (Access = private)
       xv
       yv
    end
   
    methods (Access = public)
        
        function obj = StructuredMesh(cParams)
            obj.init(cParams);
            s.coord  = obj.createCoordinates();
            s.connec = obj.createConnectivities();
            obj.createMesh(s);
        end

        function [xL,cells] = obtainLocalFromGlobalCoord(obj,xG)
            s.mesh     = obj;
            s.points.x = xG(:,1);
            s.points.y = xG(:,2);
            cellFinder = CellFinderInStructuredMesh(s);
            xL    = cellFinder.naturalCoord;
            cells = cellFinder.cells;            
        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            nDim = length(cParams.coords);
            obj.divs = cell(nDim,1);
            [obj.divs{:}] = ndgrid(cParams.coords{:});

            obj.nDivs = zeros(nDim,1);
            for i=1:nDim
                obj.nDivs(i) = size(obj.divs{1},i);
            end
        end

        function coord = createCoordinates(obj)
            nDim = length(obj.divs);
            nNodes = prod(obj.nDivs,"All");
            coord = zeros(nNodes,nDim);
            for i=1:nDim
                coord(:,i) = reshape(pagetranspose(obj.divs{i}),1,[]);
            end
        end
        
        function connec = createConnectivities(obj)
            Nx = obj.nx;
            Ny = obj.ny;
            lsw = obj.createVertex(1:Nx-1,1:Ny-1);
            lse = obj.createVertex(2:Nx,1:Ny-1);
            lne = obj.createVertex(2:Nx,2:Ny);
            lnw = obj.createVertex(1:Nx-1,2:Ny);
            connec = [lsw lse lne lnw];
        end
        
        function v = createVertex(obj,xv,yv)
            Ny = obj.ny;
            v = bsxfun(@(x,y) x + Ny*(y-1),xv',yv);
            v = v(:);
        end

        function createMesh(obj,s)
            s.kFace = 0;
            obj.mesh = Mesh.create(s);
        end
        
    end
    
end