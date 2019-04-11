classdef Interpolator < handle
    
    properties (Access = private)
       mesh
       interpolation
       grid
       zGrid
       zInterp
       x
       y
       npoints
       naturalCoord
       cellOwner
    end
    
    
    methods (Access = public)
    
        function obj = Interpolator(cParams)
            obj.init(cParams);
            obj.createMesh();
            obj.createInterpolation()
        end
        
        function z = interpolate(obj,x,y)
            obj.x = x;
            obj.y = y;
            obj.npoints = size(x,1);
            obj.obtainNaturalCoordinates();
            obj.obtainInterpolationValues();
            z = obj.zInterp;
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.grid.x = cParams.xG;
            obj.grid.y = cParams.yG;
            obj.zGrid  = cParams.zG;                 
        end
        
        function coord = createCoordinates(obj)
           coord(:,1) = reshape(obj.grid.x',1,[]);
           coord(:,2) = reshape(obj.grid.y',1,[]);
        end
        
        function createMesh(obj)
            coord  = obj.createCoordinates();
            connec = obj.createConnectivities();
            obj.mesh = Mesh();
            obj.mesh.create(coord,connec);
        end
        
        function createInterpolation(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            obj.interpolation = int;    
        end
        
        function obtainNaturalCoordinates(obj)
            xi(:,1) = obj.grid.x(1,:);
            yi(:,1) = obj.grid.y(:,1);
            [txi,xLeft] = obj.obtainCoordinate(obj.x,xi);
            [eta,yLeft] = obj.obtainCoordinate(obj.y,yi);
            obj.naturalCoord(1,:) = txi;
            obj.naturalCoord(2,:) = eta;
            nx = length(xi);
            obj.cellOwner = (xLeft) + yLeft*(nx-1);
        end
        
        function [txi,leftIndex] = obtainCoordinate(obj,xpoints,xi)
            nx = length(xi);
            dif = bsxfun(@(x, xi) (x-xi),xpoints,xi')';
            [~,imin] = min(abs(dif));
            idx = imin + nx*[0:obj.npoints-1];
            minPos = dif(idx) >= 0;
            minNeg = dif(idx) < 0;
            
            leftIndex = zeros(size(imin));  
            leftIndex(minPos) = imin(minPos);
            leftIndex(minNeg) = imin(minNeg)-1;
            
            rightIndex = zeros(size(imin));  
            rightIndex(minPos) = imin(minPos) +1;
            rightIndex(minNeg) = imin(minNeg);
            
            inc = xpoints - xi(leftIndex);
            h = xi(rightIndex) - xi(leftIndex);
            txi = 2*inc./h-1;
        end
        
        function connec = createConnectivities(obj)
            nx = length(obj.grid.x);
            ny = length(obj.grid.y);
            sw = obj.createVertex(1:nx-1,1:ny-1);
            se = obj.createVertex(2:nx,1:ny-1);
            ne = obj.createVertex(2:nx,2:ny);
            nw = obj.createVertex(1:nx-1,2:ny);
            connec = [sw se ne nw];
        end
        
        function v = createVertex(obj,xv,yv)
            ny = length(obj.grid.y);
            v = bsxfun(@(x,y) x + ny*(y-1),xv',yv);
            v = v(:);
        end
        
        function obtainInterpolationValues(obj)
            nC = obj.naturalCoord;
            obj.interpolation.computeShapeDeriv(nC);
            shapes = obj.interpolation.shape;
            z = obj.zGrid(:);
            v = zeros(obj.npoints,1);
            for inode = 1:size(shapes,1)
                nodes = obj.mesh.connec(obj.cellOwner,inode);
                znode = z(nodes);
                v(:,1) = v(:,1) + znode.*shapes(inode,:)';
            end
            obj.zInterp = v;
        end
        
    end
    
    
end