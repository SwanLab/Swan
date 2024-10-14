classdef P1Refiner < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fp1D
        fP1
    end
    
    methods (Access = public)
        
        function obj = P1Refiner(fP1)
            dofConnec = obj.computeDofConnec(fP1);
            dofCoord  = obj.computeAllDiscValues(fP1.mesh.coord,fP1.mesh.coord,fP1.mesh.connec,fP1.mesh);
            ndimf = fP1.ndimf;
            obj.fP1 = fP1;
            obj.fp1D = P1DiscontinuousFunction.create(fP1.mesh,dofConnec,dofCoord,ndimf);            
        end

        function  fP1D = compute(obj)
            fP1D = obj.fp1D;
            fP1 = obj.fP1;
            fP1D.fValues = obj.computeAllValues(fP1.fValues,fP1.getDofCoord(),fP1.mesh.connec,fP1.mesh); %%% HEre!
        end

    end
    
    methods (Access = private)

        function interpolateValues(obj)

        end
        
        function connec = computeDofConnec(obj,f)
            oldDofs = f.getDofConnec();
            newDofs = obj.computeNewDofs(f);

            vertexInCell = oldDofs;
          
            ndimf = f.ndimf;
            for iNode = 1:3
                iDof  = (iNode-1)*ndimf+(1:ndimf);
                nV(:,:,iNode) = vertexInCell(:,iDof);
            end
            nV1 = nV(:,:,1)';
            nV2 = nV(:,:,2)';
            nV3 = nV(:,:,3)';

            edgeInCell1 = squeeze(newDofs(:,1,:));
            edgeInCell2 = squeeze(newDofs(:,2,:));
            edgeInCell3 = squeeze(newDofs(:,3,:));
    
            e1d1 = edgeInCell1(1,:);
            e1d2 = edgeInCell1(2,:);
            e1d3 = edgeInCell1(3,:);

            e2d1 = edgeInCell2(1,:);
            e2d2 = edgeInCell2(2,:);
            e2d3 = edgeInCell2(3,:);

            e3d1 = edgeInCell3(1,:);
            e3d2 = edgeInCell3(2,:);
            e3d3 = edgeInCell3(3,:);

            connec(:,1,:) = [nV1 ; e1d1 ;e2d3];
            connec(:,2,:) = [e1d3; nV2 ;e3d1];
            connec(:,3,:) = [e1d2; e3d2; e2d2];
            connec(:,4,:) = [e2d1; e3d3; nV3];

            connec = reshape(connec,size(connec,1),[])';

        end

        function fEdges = computeFinEdges(obj,fDisc,connec,coordD,mesh)
            nodesDisc = obj.createDiscConnec(connec);
            
            s.nodesByElem = reshape(nodesDisc',[],mesh.nelem)';
            s.type = mesh.type;
            edge = EdgesConnectivitiesComputer(s);
            edge.compute();
                        
            
            
            s.coord  = coordD;
            s.connec = edge.nodesInEdges;
            s.kFace  = mesh.kFace -1;
            eM = Mesh.create(s); 


            s.edgeMesh = eM;
            s.fNodes   = fDisc;
            eF         = EdgeFunctionInterpolator(s);
            fInEdges = eF.compute()';           

            for iDim = 1:size(fInEdges,2)
                xc = fInEdges(:,iDim);
                xc = repmat(xc',3,1); 
                fEdges(:,iDim) = reshape(xc,[],1);
            end

        end

        function fDisc = computeDiscontinousFunction(obj,fCont,connec)
            nodesCont = reshape(connec',1,[]);
            nodesDisc = obj.createDiscConnec(connec);
            fDisc(nodesDisc,:) = fCont(nodesCont,:);
        end
        
        function nodesDisc = createDiscConnec(obj,connec)
            nNode  = size(connec,2);
            nElem  = size(connec,1);
            nodesDisc = 1:nNode*nElem;            
         %   nDime  = size(fCont,2);
            

        end


        function allF = computeAllDiscValues(obj,fCont,coord,connec,mesh)
            fDisc  = obj.computeDiscontinousFunction(fCont,connec);
            allF = obj.computeAllValues(fDisc,coord,connec,mesh);
        end

        function allFvalues = computeAllValues(obj,fDisc,coord,connec,mesh)
            coordD = obj.computeDiscontinousFunction(coord,connec);                        
            fEdges = obj.computeFinEdges(fDisc,connec,coordD,mesh);
            allFvalues = [fDisc;fEdges];
        end

            

        function newDofs = computeNewDofs(obj,f)
            dofsF = f.getDofConnec();
            maxDof = max(dofsF(:));
            nElem  = f.mesh.nelem;
            nEdges = f.mesh.edges.nEdgeByElem;
            newDofsInEdge = 3*f.ndimf;
            nNewDofs = nElem*nEdges*newDofsInEdge;
            newDofs  = maxDof + (1:nNewDofs);
            newDofs  = reshape(newDofs,nEdges,newDofsInEdge,nElem);
        end
        
    end
    
end