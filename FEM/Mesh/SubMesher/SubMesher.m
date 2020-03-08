classdef SubMesher < handle
    
   properties (Access = private)
      mesh 
      connecLocal
      eConnec
      newConnec  
      newCoord
      subMesh
   end
    
   methods (Access = public)
       
       function obj = SubMesher()
       end
        
       function s = computeSubMesh(obj,mesh)
          obj.init(mesh) 
          obj.computeNewNodes(); 
          obj.computeExtendedConnectivities();
          obj.computeConnecAndGlobalConnec();
          obj.computeNewCoordinates();
          obj.createSubMesh();
          s = obj.subMesh;
       end
    
   end
   
   methods (Access = private)

       function init(obj,mesh)
           obj.mesh = mesh;
           obj.connecLocal = [1 2 5;
                              2 3 5;
                              3 4 5;
                              1 5 4];
       end
       
       function computeExtendedConnectivities(obj)
           oldConnec = obj.mesh.connec;
           nnode     = obj.mesh.nnode;
           nnodeE = nnode + 1;
           newNodes = obj.computeNewNodes();
           connecE(:,1:nnodeE) = [oldConnec newNodes];
           obj.eConnec = connecE;
       end
       
       function nNodes = computeNewNodes(obj)
           nelemOld = obj.mesh.nelem;
           lastNode = max(obj.mesh.connec(:));
           nNodes(:,1) = (lastNode+1):(lastNode + nelemOld);
       end       
       
       function computeConnecAndGlobalConnec(obj) 
           nnodeN   = size(obj.connecLocal,2);
           nsubCells = size(obj.connecLocal,1);           
           nelemO = obj.mesh.nelem;  
           nelemN = nsubCells*nelemO;
           nConnec = zeros(nelemN,nnodeN);
           for isubcell = 1:nsubCells
              element = obj.computeGlobalElements(isubcell,nsubCells,nelemN);
              for inode = 1:nnodeN
                 node = obj.connecLocal(isubcell,inode);
                 nConnec(element,inode) = obj.eConnec(:,node);
              end
           end   
           obj.newConnec = nConnec;
       end
       
       function computeNewCoordinates(obj)
           xC = obj.computeBaricenterElementCoord();
           nCoord = [obj.mesh.coord; xC' ];    
           obj.newCoord = nCoord;
       end
       
       function xC = computeBaricenterElementCoord(obj)
           m = obj.mesh;
           q = Quadrature.set(m.geometryType);
           q.computeQuadrature('CONSTANT');
           xCenter = m.computeXgauss(q);
           xC = squeeze(xCenter);
       end
       
       function createSubMesh(obj)
           s.connec = obj.newConnec;
           s.coord  = obj.newCoord;
           m = Mesh().create(s);
           obj.subMesh = m;
       end
       
   end
    
   methods (Access = private, Static)
       
       function gElement = computeGlobalElements(isub,nsubCells,nelem)           
           gElement = isub:nsubCells:nelem ;         
       end
       
   end
    
end