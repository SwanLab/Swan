classdef MeshGenerator < handle
    
    properties (Access = public)
        mesh
    end
    
    properties (Access = private)
        coord
        connec
        type
        meshType
    end
    
    properties (Access = private)
        nElem
        columnLength
    end
    
    methods (Access = public)
        
        function obj = MeshGenerator(cParams)
            obj.init(cParams)
            obj.createCoordinates();
            obj.createConnectivity();
            obj.createMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem = cParams.nElem;
            obj.columnLength = cParams.columnLength;
            obj.type = cParams.type;
            obj.meshType = cParams.meshType;
        end
        
        function createCoordinates(obj)
            nnod = obj.nElem + 1;
            switch obj.meshType
                case 'Structured'
                    obj.coord = linspace(0,obj.columnLength,nnod)';
                case 'Unstructured'
                    x = [0;rand(nnod-2,1);1]*obj.columnLength;
                    x = sort(x);
                    obj.coord = x;
                otherwise
                    error('Invalid Mesh Type.')
            end    
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
        
        function createMesh(obj)
            s.coord  = obj.coord;  
            s.connec = obj.connec;
            s.type = obj.type;
            m = Mesh(s);
            obj.mesh = m;
        end
        
    end
    
end