classdef DimensionVariables < handle
    
    properties (Access = public)
        nnode
        nunkn
        nstre
        ndof
        nelem
        ndofPerElement 
        ngaus
        nentries
        ndim
        nunknPerField
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
            obj.nunkn          = obj.createNUnknPerField();
            obj.nstre          = obj.createNstre();
            obj.ndof           = obj.mesh.npnod*obj.nunkn;
            obj.nelem          = obj.mesh.nelem;
            obj.ndofPerElement = obj.nnode*obj.nunkn;
            obj.nentries       = obj.nelem*(obj.ndofPerElement)^2;
            obj.ndim           = obj.createNdim();
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
        
        function nUnkn = createNUnknPerField(obj) % createNUnknPerField
            switch obj.pdim
                case '2D'
                    nUnkn = 2;
                case '3D'
                    nUnkn = 3;
                case 'FILTER'
                    nUnkn = 1;
            end
        end

        function ndim = createNdim(obj)
            switch obj.pdim
                case '2D'
                    ndim = 2;
                case '3D'
                    ndim = 3;
                case 'FILTER'
                    ndim = [];
            end
        end

        function nstre = createNstre(obj)
            switch obj.pdim
                case '2D'
                    nstre = 3;
                case '3D'
                    nstre = 6;
                case 'FILTER'
                    nstre = [];
            end
        end

    end
    
end
