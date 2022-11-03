classdef DataCreator < handle
    
    properties (Access = public)
        mesh
        dim
        freeNodes
    end

    properties (Access = private)
        coord
        connec
    end

    properties (Access = private)
        nElem
        columnLength
    end

    methods
        function obj = DataCreator(cParams)
            obj.init(cParams);
            obj.createCoordinates();
            obj.createConnectivity();
            obj.createMesh();
            obj.createDimensions();
            obj.createBoundaryConditions();
        end   
    end

    methods
        function init(obj,cParams)
            obj.nElem        = cParams.nElem;
            obj.columnLength = cParams.columnLength;
        end

        function createMesh(obj)
            s.coord  = obj.coord;  
            s.connec = obj.connec;
            s.type = 'LINE';
            m = Mesh(s);
            obj.mesh = m;
        end

        function createCoordinates(obj)
            nnod = obj.nElem + 1;
%              x = [0;rand(nnod-2,1);1]*obj.columnLength;
%              x = sort(x);
%               coord = x; 
             obj.coord = linspace(0,obj.columnLength,nnod)';
        end

        function createConnectivity(obj)
            nNode = 2;
            Tnod = zeros(obj.nElem,nNode);
            e = 1;
            for iElem = 1: obj.nElem
                Tnod(iElem,1) = e;
                e = e + 1;
                Tnod(iElem,2) = e;
            end
            obj.connec = Tnod;
        end

        function createDimensions(obj)
            s.mesh = obj.mesh;
            s.pdim = '2D';
            s.ngaus = 2;
            s.type = 'Vector';
            s.name = 'x';
            s.ndimf = 2;
            s.fieldName = 'u';
            d = DimensionVariables.create(s);
            d.compute();
            obj.dim = d;
        end

        function createBoundaryConditions(obj)
            d = obj.dim;
            fixnodes = union([1,2], [d.ndofs-1,d.ndofs]);
            nodes = 1:d.ndofs;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
        end
    end
end