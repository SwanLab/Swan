classdef DimensionVector < handle
    
    properties (Access = public)
        ndofs
        nnodes
    end

    properties (GetAccess = public, SetAccess = private)
        fieldName
        scalarFields
        ndimf
        nnodeElem
        ndofsElem
    end
    
    properties (Access = private)
        mesh
    end
    
    methods (Access = public)

        function obj = DimensionVector(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            for i = 1:obj.ndimf
                name = append('field', int2str(i));
                s.mesh = obj.mesh;
                s.name = name;
                obj.scalarFields.(name) = DimensionScalar(s);
            end
            obj.ndofs     = obj.ndimf*obj.mesh.nnodes;
            obj.nnodes    = obj.mesh.nnodes;
            obj.nnodeElem = obj.mesh.interpolation.nnode;
            obj.ndofsElem = obj.nnodeElem*obj.ndimf;
        end

        function createFromScalars(obj, dims)
            obj.ndimf = numel(dims);
            for i = 1:obj.ndimf
                dim = dims{i};
                name = dim.name;
                obj.scalarFields.(name) = dim;
            end
        end
        
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh      = cParams.mesh;
            obj.ndimf = cParams.ndimf;
%             obj.fieldName = cParams.fieldName;
%             obj.interpolation = cParams.interpolation;
        end

    end

end

