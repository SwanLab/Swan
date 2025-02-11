classdef CannonicalMeshComputer < handle
    
    properties (Access = private)
        newNodes
        newCoord
        newConnec
    end
    
    properties (Access = private)
        remainingNodes
        mesh
    end
    
    methods (Access = public)
        
        function obj = CannonicalMeshComputer(cParams)
            obj.init(cParams)
        end

        function m = compute(obj)
            obj.computeNewNodes();
            obj.computeNewCoord();
            obj.computeNewConnec();
            m = obj.createCanonicalMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.remainingNodes = cParams.remainingNodes;
        end

        function computeNewNodes(obj)
            obj.newNodes = 1:length(obj.remainingNodes);
        end

        function computeNewCoord(obj)
            coord  = obj.mesh.coord;
            rNodes = obj.remainingNodes;
            nodes  = obj.newNodes;
            nCoord(nodes,:) = coord(rNodes,:);
            obj.newCoord = nCoord;
        end

        function computeNewConnec(obj)
            s.oldNodes = obj.remainingNodes;
            s.newNodes = obj.newNodes;
            c = ConnecRenumbering(s);
            connec = c.renumber(obj.mesh.connec);
            obj.newConnec = connec;
        end
        
        function m = createCanonicalMesh(obj)
            s.connec = obj.newConnec;
            s.coord  = obj.newCoord;
            s.kFace  = obj.mesh.kFace;
            m = Mesh.create(s);
        end
        
    end
    
end