classdef DimensionVector < handle
    
    properties (GetAccess = public, SetAccess = private)
        scalarFields
        nelem
        npnod
        nnode
        ndimField
        ndofPerElement
        ndof
        nstre

        ngaus % should NOT be here
    end
    
    properties (Access = private)
        mesh
    end
    
    methods (Access = public)

        function obj = DimensionVector(cParams)
            obj.init(cParams);
        end

        function create(obj, cParams)
            ndimf = cParams.ndimf;
            msh   = cParams.mesh;
            fieldName = cParams.fieldName;
            for i = 1:ndimf
                name = append(fieldName, int2str(i));
                s.mesh = msh;
                s.name = name;
                obj.scalarFields.(name) = DimensionScalar(s);
            end
            % hmmm, its the same as dimensionScalar
            obj.mesh           = msh;
            obj.nelem          = msh.nelem;
            obj.npnod          = msh.npnod;
            obj.nnode          = msh.interpolation.nnode;
            obj.ndimField      = ndimf;
            obj.ndofPerElement = obj.nnode * obj.ndimField;
            obj.ndof           = ndimf*obj.npnod;
            obj.nstre          = obj.computeNstress(); %nvoigt?
        end

        function createFromScalars(obj, dims)
            obj.ndimField = numel(dims);
            for i = 1:obj.ndimField
                dim = dims{i};
                name = dim.name;
                obj.scalarFields.(name) = dim;
            end
        end

        function applyNgaus(obj, ngaus)
            obj.ngaus = ngaus;
        end
        
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh = cParams.mesh;
        end

        function nstre = computeNstress(obj) % nvoigt?
            switch obj.ndimField
                case 1
                    nstre = 2; % ?
                case 2
                    nstre = 3;
                case 3
                    nstre = 6;
            end
        end

    end

end

