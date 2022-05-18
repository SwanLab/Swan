classdef FieldDimensions < handle

    properties (Access = public)
        ndimf
        nnodes
        ndofs
        nnodeElem
        ndofsElem
    end

    properties (Access = private)
        coord
        interpolation
    end

    methods (Access = public)

        function obj = FieldDimensions(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.nnodes    = size(obj.coord, 1);
            obj.ndofs     = obj.nnodes*obj.ndimf;
            obj.nnodeElem = obj.interpolation.nnode;
            obj.ndofsElem = obj.nnodeElem*obj.ndimf;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.coord         = cParams.coord;
            obj.ndimf         = cParams.ndimf;
            obj.interpolation = cParams.interpolation;
        end

    end

end