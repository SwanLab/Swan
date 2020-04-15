classdef SubMesher < handle
    
    properties (Access = private)
        connecLocal
        coordLocal
        coordLocalByElem        
        eConnec
        newConnec
        newCoord     
        globalToLocal
    end
    
    properties (Access = private)
        mesh
        lastNode
    end
    
    methods (Access = public)
        
        function obj = SubMesher(cParams)
            obj.init(cParams);
            obj.computeCoordLocalByElem();            
            obj.computeNewNodes();
            obj.computeExtendedConnectivities();
            obj.computeConnecAndGlobalConnec();
            obj.computeNewCoordinates();
        end
        
        function m = computeSubMesh(obj)
            s.connec = obj.newConnec;
            s.coord  = obj.newCoord;
            m = Mesh().create(s);
        end
        
        function fC = projectToSubMesh(obj,f)
            fElem = obj.computeFnodesByelem(f);
            fC    = obj.computeFinCenterElem(fElem);
        end
        
        function xE = projectToIsoSubMesh(obj,cutCells,xIsoCutCells,fullCells)

            nFull = size(fullCells,1);
            nCut  = size(cutCells,1);
            nElem = nFull + nCut;
            iFull = 1:nFull;
            iCut  = nFull + (1:nCut);             
            
            allSubCells = zeros(nElem,1);
            
            allSubCells(iFull,1) = fullCells;
            allSubCells(iCut,1)  = cutCells;
            
            xNodalAllIso = obj.computeXnodesSubCell(allSubCells);
            
            nDim  = size(xIsoCutCells,1);
            nNode = size(xIsoCutCells,2);
            xIsoAll = zeros(nDim,nNode,nElem);
            xIsoFull = xNodalAllIso(:,:,iFull);
            xIsoAll(:,:,iFull) = xIsoFull;
            xIsoAll(:,:,iCut)  = xIsoCutCells;
            
            s.connec = obj.connecLocal;
            s.coord  = obj.coordLocal;
            m = Mesh().create(s);
            int = Interpolation.create(m,'LINEAR');                        
            xE = int.interpolateFunction(xIsoAll,xNodalAllIso);
        end
        
        function xIso = computeXnodesSubCell(obj,localSubCells)            
            localSubCells = obj.globalToLocal(localSubCells);           
            xIso = obj.coordLocalByElem(:,:,localSubCells);                         
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
            obj.computeGlobalToLocal();
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
            nelemB = obj.mesh.nelem;
            nelemN = nsubCells*nelemB;
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
        
        function fC = computeFinCenterElem(obj,fElem)
            m = obj.mesh;
            q = Quadrature.set(m.geometryType);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(xV);
            shapes = int.shape;
            fCenter = m.interpolateFunction(shapes,fElem);
            fC = squeeze(fCenter);
        end
        
        function fElem = computeFnodesByelem(obj,f)
            nNode  = obj.mesh.nnode;
            nDime  = size(f,2);
            nElem  = obj.mesh.nelem;
            coordE = zeros(nDime,nNode,nElem);
            coords = f';
            for inode = 1:nNode
                nodes = obj.mesh.connec(:,inode);
                coordNodes = coords(:,nodes);
                coordE(:,inode,:) = coordNodes;
            end
            fElem = coordE;
        end
        
        function subCellIso = computeSubCellIso(obj,cellSubMesh)
            
            subCellIso = cell(cellSubMesh);
        end
        
        function computeGlobalToLocal(obj)
            bConnec = obj.mesh.connec;
            nnode   = size(bConnec,2);
            nElem   = size(bConnec,1);
            cell = repmat((1:nnode)',1,nElem);
            obj.globalToLocal = cell(:);            
        end
        
        function computeCoordLocalByElem(obj)
            nodes = obj.connecLocal(:,:);
            coord = obj.coordLocal;
            nNode = size(nodes,2);
            nElem = size(nodes,1);
            nDime = size(coord,2);
            coords = transpose(coord);
            x = zeros(nDime,nNode,nElem);
            for inode = 1:nNode
                node = nodes(:,inode);
                x(:,inode,:) = coords(:,node);
            end
            obj.coordLocalByElem = x;
        end
        
    end
    
    methods (Access = private, Static)
        
        function gElement = computeGlobalElements(isub,nsubCells,nelem)
            gElement = isub:nsubCells:nelem ;
        end
        
    end
    
end