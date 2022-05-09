classdef DimensionVector < handle
    
    properties (GetAccess = public, SetAccess = private)
        fieldName
        scalarFields
        nnodes
        nnodeElem
        ndimField
        ndofPerElement
        ndof
    end
    
    properties (Access = private)
        mesh
    end
    
    methods (Access = public)

        function obj = DimensionVector(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            ndimf = obj.ndimField;
            for i = 1:ndimf
                name = append(obj.fieldName, int2str(i));
                s.mesh = obj.mesh;
                s.name = name;
                obj.scalarFields.(name) = DimensionScalar(s);
            end
            obj.ndof           = ndimf*obj.mesh.nnodes;
            obj.nnodes         = obj.mesh.nnodes;
            obj.nnodeElem      = obj.mesh.interpolation.nnode;
            obj.ndofPerElement = obj.nnodeElem*obj.ndimField;
        end

        function createFromScalars(obj, dims)
            obj.ndimField = numel(dims);
            for i = 1:obj.ndimField
                dim = dims{i};
                name = dim.name;
                obj.scalarFields.(name) = dim;
            end
        end
        
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh      = cParams.mesh;
            obj.ndimField = cParams.ndimf;
            obj.fieldName = cParams.fieldName;
        end

    end

end

