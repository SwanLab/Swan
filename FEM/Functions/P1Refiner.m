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
            dofCoord  = obj.computeDiscFunction(fP1.mesh.coord,fP1.mesh.connec,fP1.mesh);
            ndimf = fP1.ndimf;
            obj.fP1 = fP1;
            obj.fp1D = P1DiscontinuousFunction.create(fP1.mesh,dofConnec,dofCoord,ndimf);            
        end

        function  fP1D = compute(obj)
            fP1D = obj.fp1D;
            fP1 = obj.fP1;
            fP1D.fValues = obj.computeDiscFunction(fP1.fValues,fP1.mesh.connec,fP1.mesh); %%% HEre!
        end

    end
    
    methods (Access = private)

        function interpolateValues(obj)

        end
        
        function connec = computeDofConnec(obj,f)
            oldDofs = f.getDofConnec();
            newDofs = obj.computeNewDofs(f);

            vertexInCell = oldDofs;
            nV1 = vertexInCell(:,1)';
            nV2 = vertexInCell(:,2)';
            nV3 = vertexInCell(:,3)';

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

        function newCoord = createNewValues(obj,fDisc,nodesDisc,mesh)
            
            s.nodesByElem = reshape(nodesDisc',[],mesh.nelem)';
            s.type = mesh.type;
            edge = EdgesConnectivitiesComputer(s);
            edge.compute();
                        
            
            
            s.coord  = fDisc;
            s.connec = edge.nodesInEdges;
            s.kFace  = mesh.kFace -1;
            eM = Mesh.create(s); 


            s.edgeMesh = eM;
            s.fNodes   = fDisc;
            eF         = EdgeFunctionInterpolator(s);
            newEdgeCoord = eF.compute()';           

            for iCoord = 1:size(newEdgeCoord,2)
                xc = newEdgeCoord(:,iCoord);
                xc = repmat(xc',3,1); 
                newCoord(:,iCoord) = reshape(xc,[],1);
            end

        end

        function allCoord = computeDiscFunction(obj,fCont,connec,mesh)

            nNode  = size(connec,2);
            nElem  = size(connec,1);
         %   nDime  = size(fCont,2);
            nodesCont = reshape(connec',1,[]);
            nodesDisc = 1:nNode*nElem;
            fDisc(nodesDisc,:) = fCont(nodesCont,:);
        %    fVals = reshape(fDisc',nDime,nNode,[]);


            newCoord = obj.createNewValues(fDisc,nodesDisc,mesh);

            allCoord = [fDisc;newCoord];
        end

            

        function newDofs = computeNewDofs(obj,f)
            dofsF = f.getDofConnec();
            maxDof = max(dofsF(:));
            nElem  = f.mesh.nelem;
            nEdges = f.mesh.edges.nEdgeByElem;
            newDofsInEdge = 3;
            nNewDofs = nElem*nEdges*newDofsInEdge;
            newDofs  = maxDof + (1:nNewDofs);
            newDofs  = reshape(newDofs,nEdges,newDofsInEdge,nElem);
        end
        
    end
    
end