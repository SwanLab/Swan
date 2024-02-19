classdef Geometry < handle
    
    properties (GetAccess = public, SetAccess = protected)
        dvolu
    end

    properties (SetAccess = private, GetAccess = protected)
        xFE
    end

    methods (Access = public, Static)

        function obj = create(cParams)
            f = GeometryFactory();
            obj = f.create(cParams);
        end

    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.xFE = cParams.xFE;
        end

    end

end