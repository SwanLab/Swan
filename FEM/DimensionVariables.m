classdef DimensionVariables < handle
    
    properties (Access = public)
        nnodeElem
        ndimField
        nstre
        ndof
        ndofPerElement
        ndim
        nnodes
    end

    properties (Access = private)
        pdim
        mesh
    end
    
    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.type
                case 'Scalar'
                    obj = DimensionScalar(cParams);
                case 'Vector'
                    obj = DimensionVector(cParams);
            end
        end
    end
    
    methods (Access = public)

        function obj = DimensionVariables(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.nnodeElem      = obj.mesh.nnodeElem;
            obj.ndimField      = obj.createDimPerField();
            obj.ndof           = obj.mesh.nnodes*obj.ndimField;
            obj.ndofPerElement = obj.nnodeElem*obj.ndimField;
            obj.ndim           = obj.createNdim();
            obj.nnodes          = obj.mesh.nnodes;
        end

        function applyNdimfield(obj, num)
            obj.ndimField = num;
            obj.ndofPerElement = obj.nnodeElem*obj.ndimField;
            obj.ndof           = obj.mesh.nnodes*obj.ndimField;
        end

    end
    
    methods (Access = private)

        function obj = init(obj, cParams)
            obj.mesh  = cParams.mesh;
            obj.pdim  = cParams.pdim;
        end
        
        function ndimf = createDimPerField (obj) % createNUnknPerField
            switch obj.pdim
                case '1D'
                    ndimf = 1;
                case '2D'
                    ndimf = 2;
                case '3D'
                    ndimf = 3;
            end
        end

        function ndim = createNdim(obj)
            switch obj.pdim
                case '1D'
                    ndim = 1;
                case '2D'
                    ndim = 2;
                case '3D'
                    ndim = 3;
            end
        end

        function nstre = createNstre(obj)
            switch obj.pdim
                case '1D'
                    nstre = 2;
                case '2D'
                    nstre = 3;
                case '3D'
                    nstre = 6;
            end
        end

    end
    
end
