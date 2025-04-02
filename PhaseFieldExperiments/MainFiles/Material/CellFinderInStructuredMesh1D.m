classdef CellFinderInStructuredMesh1D < handle
    
    properties (Access = public)
        naturalCoord
        cells
    end

    properties (Access = private)
        sMesh
        xLeftIndex
        points
    end

    methods (Access = public)

        function obj = CellFinderInStructuredMesh1D(cParams)
            obj.init(cParams);
            obj.obtainLeftIndeces();
            obj.obtainCellNumber();
            obj.obtainNaturalCoordinates();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.sMesh  = cParams.mesh;
            obj.points = cParams.points;
        end

        function obtainLeftIndeces(obj)
            xi = obj.sMesh.x;
            obj.xLeftIndex = obj.obtainLeftIndex(obj.points.x,xi);
        end

        function obtainCellNumber(obj)
            obj.cells(:,1) = obj.xLeftIndex;        
        end

        function obtainNaturalCoordinates(obj)
            xi = obj.sMesh.x;
            xL = obj.xLeftIndex;
            txi = obj.obtainCoordinate(obj.points.x,xi,xL);

            obj.naturalCoord = txi;
        end

    end

    methods (Access = private, Static)

        function lIndex = obtainLeftIndex(xpoints,xi)
            nx = length(xi);
            npoints = length(xpoints);
            dif = bsxfun(@(x, xi) (x-xi),xpoints,xi')';
            [~,imin] = min(abs(dif));
            idx = imin + nx*(0:npoints-1);
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
