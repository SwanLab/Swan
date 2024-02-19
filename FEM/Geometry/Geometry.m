classdef Geometry < handle

    properties (SetAccess = private, GetAccess = protected)
        xFE
        coord
        interpolation
    end

    methods (Access = public, Static)

        function obj = create(cParams)
            f = GeometryFactory();
            obj = f.create(cParams);
        end

    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.coord = cParams.coord;
            obj.xFE = cParams.xFE;
        end

    end

end