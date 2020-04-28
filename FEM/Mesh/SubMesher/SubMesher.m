classdef SubMesher < handle
    
    properties (Access = private)
        connecLocal
        coordLocal
        newConnec
        newCoord
    end
    
    properties (Access = private)
        mesh
        lastNode
    end
    
    methods (Access = public)
        
        function obj = SubMesher(cParams)
            obj.init(cParams);
            obj.computeNewConnec();
            obj.computeNewCoordinates();
        end
        
        function m = computeSubMesh(obj)
            s.connec = obj.newConnec;
            s.coord  = obj.newCoord;
            m = Mesh().create(s);
        end
        
        function m = computeLocalSubMesh(obj,localSubCells)
            s.coord  = obj.coordLocal;
            s.connec = obj.connecLocal(localSubCells,:);
            m = Mesh().create(s);            
        end
              
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.lastNode    = cParams.lastNode;
            obj.connecLocal = [1 2 5;
                2 3 5;
                3 4 5;
                1 5 4];
            obj.coordLocal = obj.computeCoordLocal();
        end
        
        function xIsoCoord = computeCoordLocal(obj)
            xIsoQuad = obj.computeXisoQuad();
            xIsoCent = [0 0];
            xIsoCoord = [xIsoQuad;xIsoCent];
        end
        
        function x = computeXisoQuad(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            x = int.pos_nodes;
        end
        
        function computeNewConnec(obj)
            connecE = obj.computeExtendedConnectivities();
            nnodeN   = size(obj.connecLocal,2);
            nsubCells = size(obj.connecLocal,1);
            nelemB = obj.mesh.nelem;
            nelemN = nsubCells*nelemB;
            nConnec = zeros(nelemN,nnodeN);
            for isubcell = 1:nsubCells
                element = obj.computeGlobalElements(isubcell,nsubCells,nelemN);
                for inode = 1:nnodeN
                    node = obj.connecLocal(isubcell,inode);
                    nConnec(element,inode) = connecE(:,node);
                end
            end
            obj.newConnec = nConnec;
        end
        
        function connecE = computeExtendedConnectivities(obj)
            oldConnec = obj.mesh.connec;
            nnode     = obj.mesh.nnode;
            nnodeE    = nnode + 1;
            newNodes = obj.computeNewNodes();
            connecE(:,1:nnodeE) = [oldConnec newNodes];
        end
        
        function nNodes = computeNewNodes(obj)
            nelemOld = obj.mesh.nelem;
            nNodes(:,1) = (obj.lastNode+1):(obj.lastNode + nelemOld);
        end        
        
        function computeNewCoordinates(obj)
            xC = obj.mesh.computeBaricenter();
            nCoord = [obj.mesh.coord; xC' ];
            obj.newCoord = nCoord;
        end
        
    end
    
    methods (Access = private, Static)
        
        function gElement = computeGlobalElements(isub,nsubCells,nelem)
            gElement = isub:nsubCells:nelem ;
        end
        
    end
    
end