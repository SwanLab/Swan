classdef FrameSpecimenGeometry < SpecimenGeometry

    properties (Access = private)
        gapTolerance
    end

    methods (Access = public)

        function obj = FrameSpecimenGeometry(cParams)
            obj.init(cParams);
            obj.setGrips(cParams);
        end

        function ls = computeLevelSetFunction(obj,mesh)
            t = obj.hingeWidth;
            g.type        = 'RectangleInclusion';
            g.xSide       = obj.sideLength - 2*t;
            g.ySide       = obj.height - 2*t;
            g.xCoorCenter = obj.center(1);
            g.yCoorCenter = obj.center(2);
            ls = GeometricalFunction(g).computeLevelSetFunction(mesh);
        end

        function mesh = createBackgroundMesh(obj,refinement,~,~)
            h = refinement*obj.hingeWidth/4;
            W = obj.sideLength;
            H = obj.height;
            mesh = TriangleMesh(W,H,round(W/h)+1,round(H/h)+1);
        end

        function s = computeBarrierParameters(obj)
            s = computeBarrierParameters@SpecimenGeometry(obj);
            s.gapTolerance       = obj.gapTolerance;
            s.stiffnessReduction = 10;
        end

    end

    methods (Static, Access = public)

        function cP = paperFrameParams()
            cP.width         = 36;
            cP.height        = 9;
            cP.wallThickness = 2;
            cP.punchWidth    = 2;
            cP.gapTolerance  = 0.02;
        end

    end

    methods (Access = protected)

        function points = computeRefinePoints(obj)
            points = obj.shape.Vertices;
        end

        function s = computeMeshSizes(obj,refinement)
            h = refinement*obj.hingeWidth/4;
            s.sizeMin = h;
            s.sizeMax = h;
            s.distMin = obj.hingeWidth;
            s.distMax = 2*obj.hingeWidth;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            W = cParams.width;
            H = cParams.height;
            t = cParams.wallThickness;
            outer = polyshape([0 W W 0],[0 0 H H]);
            inner = polyshape([t W-t W-t t],[t t H-t H-t]);
            obj.shape          = subtract(outer,inner);
            obj.hingeWidth     = t;
            obj.gapTolerance   = cParams.gapTolerance;
            obj.cornerLigament = t;
            obj.hingeTips      = zeros(0,2);
            obj.sideLength     = W;
            obj.height         = H;
            obj.center         = [W/2,H/2];
        end

        function setGrips(obj,cParams)
            W = cParams.width;
            H = cParams.height;
            a = cParams.punchWidth/2;
            obj.bottomGrips     = [0, W, 0];
            obj.topGrips        = [W/2-a, W/2+a, H];
            obj.latticeVertices = [W/2-a, H; W/2+a, H];
        end

    end

end
