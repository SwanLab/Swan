classdef QuadToTriMeshConverter < handle
    
    properties (Access = public)
       subMesh
       localMesh
    end
    
    properties (Access = private)
        mesh
        lastNode
    end
    
    methods (Access = public)
        
        function obj = QuadToTriMeshConverter(cParams)
            obj.init(cParams);
            obj.computeLocalMesh();
            obj.computeSubMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.lastNode    = cParams.lastNode;
        end
        
        function computeLocalMesh(obj)
            s.connec = [1 2 5;
                2 3 5;
                3 4 5;
                1 5 4];
            s.coord =  obj.computeCoordLocal();
            obj.localMesh = Mesh.create(s);
        end
        
        function xIsoCoord = computeCoordLocal(obj)
            xIsoQuad = obj.computeXisoQuad();
            xIsoCent = [0 0];
            xIsoCoord = [xIsoQuad;xIsoCent];
        end
        
        function computeSubMesh(obj)
            s.connec = obj.computeSubMeshConnec();
            s.coord  = obj.computeSubMeshCoord();
            obj.subMesh = Mesh.create(s);
        end
        
        function x = computeXisoQuad(obj)
            int = Interpolation.create(obj.mesh.type,'LINEAR');
            x = int.pos_nodes;
        end
        
        function newConnec = computeSubMeshConnec(obj)
            connecE = obj.computeExtendedConnectivities();
            nnodeN   = size(obj.localMesh.connec,2);
            nsubCells = size(obj.localMesh.connec,1);
            nelemB = obj.mesh.nelem;
            nelemN = nsubCells*nelemB;
            newConnec = zeros(nelemN,nnodeN);
            for isubcell = 1:nsubCells
                element = obj.computeGlobalElements(isubcell,nsubCells,nelemN);
                for inode = 1:nnodeN
                    node = obj.localMesh.connec(isubcell,inode);
                    newConnec(element,inode) = connecE(:,node);
                end
            end
        end
        
        function connecE = computeExtendedConnectivities(obj)
            oldConnec = obj.mesh.connec;
            nnode     = obj.mesh.nnodeElem;
            nnodeE    = nnode + 1;
            newNodes = obj.computeNewNodes();
            connecE(:,1:nnodeE) = [oldConnec newNodes];
        end
        
        function nNodes = computeNewNodes(obj)
            nelemOld = obj.mesh.nelem;
            nNodes(:,1) = (obj.lastNode+1):(obj.lastNode + nelemOld);
        end
        
        function newCoord = computeSubMeshCoord(obj)
            xC = obj.mesh.computeBaricenter();
            newCoord = [obj.mesh.coord; xC' ];
        end
        
    end
    
    methods (Access = private, Static)
        
        function gElement = computeGlobalElements(isub,nsubCells,nelem)
            gElement = isub:nsubCells:nelem ;
        end
        
    end
    
end