classdef SubMesher < handle
    
   properties (Access = private)
      connecLocal
      coordLocal
      eConnec
      newConnec  
      newCoord
      subMesh
   end
   
   properties (Access = private)
      mesh 
      lastNode       
   end
    
   methods (Access = public)
       
       function obj = SubMesher(cParams)
          obj.init(cParams)            
       end
        
       function s = computeSubMesh(obj)
          obj.computeNewNodes(); 
          obj.computeExtendedConnectivities();
          obj.computeConnecAndGlobalConnec();
          obj.computeNewCoordinates();
          obj.createSubMesh();
          s = obj.subMesh;
       end
       
       function fC = projectToSubMesh(obj,f)
           fElem = obj.computeFnodesByelem(f);
           fElem   = permute(fElem,[3 2 1]);                      
           
           fC    = obj.computeFinCenterElem(fElem);
       end
       
       function xE = projectToIsoSubMesh(obj,subCellIso,xIso)
           nodes = obj.connecLocal(subCellIso,:);            
           coord = obj.coordLocal;
           nNode = size(nodes,2);
           nElem = size(nodes,1);
           nDime = size(coord,2);
           xIsoSubCell = zeros(nNode,nDime,nElem);
           for inode = 1:nNode
               node = nodes(:,inode);
               xIsoSubCell(inode,:,:) = coord(node,:)';
           end
               
           
           xIso = permute(xIso,[2 3 1]);
           s.connec = obj.connecLocal;
           s.coord  = obj.coordLocal;
           m = Mesh().create(s);
           xIsoSubCell   = permute(xIsoSubCell,[3 2 1]);                      
         
           xE = m.interpolateFunction(xIso,xIsoSubCell);
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
           nNodes(:,1) = (obj.lastNode+1):(obj.lastNode + nelemOld);
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
       
    function fC = computeFinCenterElem(obj,fElem)
           m = obj.mesh;
           q = Quadrature.set(m.geometryType);
           q.computeQuadrature('CONSTANT');
           xV = q.posgp;
           fCenter = m.interpolateFunction(xV,fElem);
           fC = squeeze(fCenter);            
       end
    
       function fElem = computeFnodesByelem(obj,f)
            nNode  = obj.mesh.nnode;
            nDime  = size(f,2);
            nElem  = obj.mesh.nelem;            
            coordE = zeros(nNode,nDime,nElem);
            coords = f';
            for inode = 1:nNode
                nodes = obj.mesh.connec(:,inode);
                coordNodes = coords(:,nodes);
                coordE(inode,:,:) = coordNodes;
            end
            fElem = coordE;               
       end       
       
   end
    
   methods (Access = private, Static)
       
       function gElement = computeGlobalElements(isub,nsubCells,nelem)           
           gElement = isub:nsubCells:nelem ;         
       end
       
   end
    
end