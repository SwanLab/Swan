classdef PolyshapeSignedDistance < handle

    properties (Access = private)
        shape
        segmentStarts
        segmentEnds
    end

    methods (Access = public)

        function obj = PolyshapeSignedDistance(shape)
            obj.shape = shape;
            obj.computeBoundarySegments();
        end

        function d = compute(obj,x)
            d = obj.computeUnsignedDistance(x);
            isInside = isinterior(obj.shape,x(:,1),x(:,2));
            d(isInside) = -d(isInside);
        end

    end

    methods (Access = private)

        function computeBoundarySegments(obj)
            p = zeros(0,2);
            q = zeros(0,2);
            for b = 1:numboundaries(obj.shape)
                [x,y] = boundary(obj.shape,b);
                p = [p;x(1:end-1),y(1:end-1)];
                q = [q;x(2:end),y(2:end)];
            end
            obj.segmentStarts = p;
            obj.segmentEnds   = q;
        end

        function d = computeUnsignedDistance(obj,x)
            d = inf(size(x,1),1);
            for i = 1:size(obj.segmentStarts,1)
                d = min(d,obj.computeDistanceToSegment(x,obj.segmentStarts(i,:),obj.segmentEnds(i,:)));
            end
        end

    end

    methods (Static, Access = private)

        function d = computeDistanceToSegment(x,p,q)
            t  = q - p;
            r  = x - p;
            xi = min(max(r*t'/(t*t'),0),1);
            d  = vecnorm(r - xi.*t,2,2);
        end

    end

end
