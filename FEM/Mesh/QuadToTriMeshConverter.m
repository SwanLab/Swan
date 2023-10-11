classdef QuadToTriMeshConverter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        localConnec
        fullCoord
        fullConnec
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function tMesh = convert(obj, qMesh)
            obj.init(qMesh);
            obj.createLocalConnectivities();
            obj.computeFullCoordinates();
            obj.computeFullConnectivities();
            tMesh = obj.createTriMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,m)
            obj.mesh = m;
        end

        function createLocalConnectivities(obj)
            conn1 = [1 2 5];
            conn2 = [2 3 5];
            conn3 = [3 4 5];
            conn4 = [4 1 5];
            obj.localConnec = [conn1;conn2;conn3;conn4];
        end

        function computeFullCoordinates(obj)
            coord    = obj.mesh.coord;
            newCoord = obj.computeNewCoordinates();
            obj.fullCoord = [coord;newCoord'];
        end

        function nCoord = computeNewCoordinates(obj)
            nCoord = obj.mesh.computeBaricenter();
        end  

        function computeFullConnectivities(obj)
            nSubElem = 4;
            nNodeTri = 3;
            nElem    = obj.mesh.nelem;
            nNodes   = obj.mesh.nnodes;
            nodes = obj.mesh.connec;
            newNodes = nNodes + (1:nElem);
            allNodes = [nodes, newNodes'];
            fConnec = zeros(nElem*nSubElem, nNodeTri);
            for iSubElem = 1:nSubElem
                localConnecI = obj.localConnec(iSubElem,:);
                index = nElem*(iSubElem-1) + (1:nElem);
                fConnec(index,:) = allNodes(:,localConnecI); 
            end
            obj.fullConnec = fConnec;
        end  
        
        function m = createTriMesh(obj)
            s.coord  = obj.fullCoord;
            s.connec = obj.fullConnec;
            m = Mesh(s);
        end
    end
    
end