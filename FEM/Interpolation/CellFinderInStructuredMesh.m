classdef CellFinderInStructuredMesh < handle
    
    
    properties (Access = public)
        naturalCoord
        cells
    end
    
    properties (Access = private)
        sMesh
        xLeftIndex
        yLeftIndex
        points
    end
    
    methods (Access = public)
        
        function obj = CellFinderInStructuredMesh(cParams)
            obj.sMesh   = cParams.mesh;
            obj.points = cParams.points;
            obj.obtainLeftIndeces();
            obj.obtainNodesOfCellsWithInterpolationPoints();
            obj.obtainNaturalCoordinates();
        end
        
    end
    
    methods (Access = private)
        
        function obtainLeftIndeces(obj)
            xi(:,1) = obj.sMesh.x(1,:);
            yi(:,1) = obj.sMesh.y(:,1);
            obj.xLeftIndex = obj.obtainLeftIndex(obj.points.x,xi);
            obj.yLeftIndex = obj.obtainLeftIndex(obj.points.y,yi);
        end        
        
        function obtainNodesOfCellsWithInterpolationPoints(obj)
            connec = obj.sMesh.mesh.connec;
            xL = obj.xLeftIndex;
            yL = obj.yLeftIndex;
            nx = obj.sMesh.nx;
            cellNumber = (xL) + (nx-1)*(yL-1);
            obj.cells = connec(cellNumber,:);
        end
        
        function obtainNaturalCoordinates(obj)
            xi(:,1) = obj.sMesh.x(1,:);
            xL = obj.xLeftIndex;
            txi = obj.obtainCoordinate(obj.points.x,xi,xL);
            
            yL = obj.yLeftIndex;
            yi(:,1) = obj.sMesh.y(:,1);
            eta = obj.obtainCoordinate(obj.points.y,yi,yL);
            
            obj.naturalCoord(1,:) = txi;
            obj.naturalCoord(2,:) = eta;
        end

    end
    
    methods (Access = private, Static)
        
        function lIndex = obtainLeftIndex(xpoints,xi)
            nx = length(xi);
            npoints = length(xpoints);
            dif = bsxfun(@(x, xi) (x-xi),xpoints,xi')';
            [~,imin] = min(abs(dif));
            idx = imin + nx*[0:npoints-1];
            minPos = dif(idx) >= 0;
            minNeg = dif(idx) < 0;            
            lIndex = zeros(size(imin));
            lIndex(minPos) = imin(minPos);
            lIndex(minNeg) = imin(minNeg)-1;
        end        
        
        function txi = obtainCoordinate(xpoints,xi,index)
            inc = xpoints - xi(index);
            h = xi(index + 1) - xi(index);
            txi = 2*inc./h-1;
        end

    end
    
end