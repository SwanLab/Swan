classdef DimensionVariables < handle
    
    properties (Access = public)
        nnode
        ndimField
        nstre
        ndof
        nelem
        ndofPerElement
        ngaus
        nentries
        ndim
        nunknPerField
        nt
        npnod
    end

    properties (Access = private)
        pdim
        mesh
    end
    
    methods (Access = public)

        function obj = DimensionVariables(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.nnode          = obj.mesh.nnode;
            obj.ndimField      = obj.createDimPerField();
            obj.nstre          = obj.createNstre();
            obj.ndof           = obj.mesh.npnod*obj.ndimField;
            obj.nelem          = obj.mesh.nelem;
            obj.ndofPerElement = obj.nnode*obj.ndimField;
            obj.nentries       = obj.nelem*(obj.ndofPerElement)^2;
            obj.ndim           = obj.createNdim();
            obj.nt             = obj.ngaus*obj.nelem*obj.nstre;
            obj.npnod          = obj.mesh.npnod;
        end

        function applyNUnknPerField(obj, num)
            obj.nunknPerField = num;
        end
        
    end
    
    methods (Access = private)

        function obj = init(obj, cParams)
            obj.mesh  = cParams.mesh;
            obj.pdim  = cParams.pdim;
            obj.ngaus = cParams.ngaus;
        end
        
        function ndimf = createDimPerField(obj) % createNUnknPerField
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
