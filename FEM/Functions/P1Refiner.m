classdef P1Refiner < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        coordD
        coordInEdges
    end
    
    properties (Access = private)
        fCoarse
        meshFine
    end
    
    methods (Access = public)
        
        function obj = P1Refiner(fCoarse,meshFine)
            obj.fCoarse = fCoarse;
            obj.meshFine = meshFine;
        end

        function  fFine = compute(obj)
            obj.coordD       = obj.computeCoordDisc();
            obj.coordInEdges = obj.computeCoordInEdges();

            s.dofCoord  = obj.computeDofCoord();
            

            s.mesh      = obj.meshFine;            
            s.dofConnec = obj.computeDofConnec();
            s.dofs.getNumberDofs = size(s.dofCoord,1);
            s.fValues = obj.computeAllFValues(); %HERE!!!
            s.order  = 'P1D';

            fFine = LagrangianFunction(s);                   
        end

    end
    
    methods (Access = private)

        function allF = computeDofCoord(obj)
          %  connec = obj.fCoarse.mesh.connec;                
            mesh      =  obj.fCoarse.mesh;            
%             nodesDisc = (obj.fCoarse.dofConnec);
            nodesDisc = obj.getDofConnecFromVector(obj.fCoarse.getDofConnec());
             
            coordD    = obj.fCoarse.getDofCoord();
            coordD = obj.coordD;%coordD;%obj.computeDiscontinousFunction(coordC,connec,nodesDisc);
            ndim = obj.fCoarse.ndimf;
      %      allF      = obj.computeAllValues(coordD,nodesDisc,mesh,ndim);
            coordEdges = obj.computeFinEdges(coordD,nodesDisc,mesh,ndim);
            allF = [obj.fCoarse.getDofCoord();coordEdges];            
        end

        function node = getDofConnecFromVector(obj,dofConnec)
          nNode = obj.fCoarse.mesh.nnodeElem;
          ndimf = obj.fCoarse.ndimf;
          for iNode = 1:nNode
            iDof   = (iNode-1)*ndimf+1;              
            node(:,iNode) = (dofConnec(:,iDof)-1)/ndimf+1;
          end
        end            

        function cD = computeCoordDisc(obj)
            cC     = obj.fCoarse.mesh.coord;
            connec = obj.fCoarse.mesh.connec;
            %nodesDisc = obj.createDiscConnec(connec); 
            nodesDisc = obj.getDofConnecFromVector(obj.fCoarse.getDofConnec());
            cD = obj.computeDiscontinousFunction(cC,connec,nodesDisc);                                
        end

        function cEdges = computeCoordInEdges(obj)
            cD     = obj.coordD;
          %  connec = obj.fCoarse.mesh.connec;  
          %  nodesDisc = obj.createDiscConnec(connec);
          %  nodesDisc = reshape(nodesDisc',[],obj.fCoarse.mesh.nelem)';

           nodesDisc = obj.getDofConnecFromVector(obj.fCoarse.getDofConnec());
            
            

            %connec = obj.fCoarse.dofConnec;
            mesh   = obj.fCoarse.mesh;            
            cEdges = obj.computeFinEdges(cD,nodesDisc,mesh,1);
        end

        function allFvalues = computeAllFValues(obj)
            fDisc     = obj.fCoarse.fValues;
           % fDisc = fDisc(:);
            nodesDisc = obj.getDofConnecFromVector(obj.fCoarse.getDofConnec);
            %nodesDisc = (obj.fCoarse.dofConnec);
            mesh      = obj.fCoarse.mesh;
      %      ndim = 1;
            fEdges = obj.computeFinEdgesWithouRep(fDisc,nodesDisc,mesh);
            for iD = 1:obj.fCoarse.ndimf
                fI = fEdges(:,iD);
                fII(:,iD,:) = reshape(fI,[],obj.fCoarse.mesh.nelem);
            end
            fEdges = fII;

            nEdge = 3;
            nElemInEdge = 3;
            for iEdge = 1:nEdge
                for iDim = 1:obj.fCoarse.ndimf
                    nod = squeeze(fII(iEdge,iDim,:));
                    idfEl = (iEdge-1)*nElemInEdge+(1:nElemInEdge);
                    nodR(idfEl,iDim,:) =repmat(nod',[nElemInEdge,1]);

                end
            end
            A = nodR; A = reshape(permute(A,[1,3,2]),[],size(A,2));
            fEdges = A;
            %fEdges = reshape(fII,[],obj.fCoarse.ndimf,obj.fCoarse.mesh.nelem);
            
            %fEdges = permute(fEdges,[2  3 1]);
            %fEdges = repmat(fEdges,[nElemInEdge 1 1]);
            %fEdges = reshape(permute(fEdges,[1,3,2]),size(fEdges,2),[]);
            %fEdges = fEdges';
            allFvalues = [obj.fCoarse.fValues;fEdges];            
%            allFvalues = obj.computeAllValues(fDisc,nodesDisc,mesh,ndim);
%            allFvalues = reshape(reshape(allFvalues,[],ndim)',[],1);
       %     allFvalues = reshape(allFvalues,[],2);
        end

        function allFvalues = computeAllValues(obj,fDisc,nodesDisc,mesh,ndim)
        %    coordD = obj.computeDiscontinousFunction(coord,connec);                        
            fEdges = obj.computeFinEdges(fDisc,nodesDisc,mesh,ndim);
            allFvalues = [fDisc;fEdges];
        end
        
        
        function dofConnec = computeDofConnec(obj)
            fCoarse = obj.fCoarse;
            oldDofs = fCoarse.getDofConnec();
            newDofs = obj.computeNewDofs(fCoarse);

            vertexInCell = oldDofs;          
            ndimf = fCoarse.ndimf;
            for iNode = 1:3
                     iDof  = (iNode-1)*ndimf+(1:ndimf);
                    nV(:,:,iNode) = vertexInCell(:,iDof);                
                
                
            end

            edgeIn = newDofs;

  
            nV1(1:ndimf,:) = nV(:,:,1)';
            %nV1 = nV1(:);
            %nV1 = nV1';


            nV2(1:ndimf,:) = nV(:,:,2)';
            %nV2 = nV2(:);
            %nV2 = nV2';

            nV3(1:ndimf,:) = nV(:,:,3)';
            %nV3 = nV3(:);
            %nV3 = nV3';


            % edgeInCell1 = squeeze(newDofs(:,1,:));
            % edgeInCell2 = squeeze(newDofs(:,2,:));
            % edgeInCell3 = squeeze(newDofs(:,3,:));            
    
            % e1d1 = edgeInCell1(1,:);
            % e1d2 = edgeInCell1(2,:);
            % e1d3 = edgeInCell1(3,:);
            % 
            % e2d1 = edgeInCell2(1,:);
            % e2d2 = edgeInCell2(2,:);
            % e2d3 = edgeInCell2(3,:);
            % 
            % e3d1 = edgeInCell3(1,:);
            % e3d2 = edgeInCell3(2,:);
            % e3d3 = edgeInCell3(3,:);

            nNewDofs = size(newDofs,4);
            e1d1(1:ndimf,:) = (squeeze(edgeIn(:,1,1,:)));
            %e1d1c = reshape(e1d1,[],nElem);            
            e1d2(1:ndimf,:) = (squeeze(edgeIn(:,2,1,:)));
            e1d3(1:ndimf,:) = (squeeze(edgeIn(:,3,1,:)));

            e2d1(1:ndimf,:) = (squeeze(edgeIn(:,1,2,:)));
            e2d2(1:ndimf,:) = (squeeze(edgeIn(:,2,2,:)));
            e2d3(1:ndimf,:) = (squeeze(edgeIn(:,3,2,:)));

            e3d1(1:ndimf,:) = (squeeze(edgeIn(:,1,3,:)));
            e3d2(1:ndimf,:) = (squeeze(edgeIn(:,2,3,:)));
            e3d3(1:ndimf,:) = (squeeze(edgeIn(:,3,3,:)));



            dofConnec(:,1,:) = [nV1 ; e1d1 ;e2d3];
            dofConnec(:,2,:) = [e1d3; nV2 ;e3d1];
            dofConnec(:,3,:) = [e1d2; e3d2; e2d2];
            dofConnec(:,4,:) = [e2d1; e3d3; nV3];

            dofConnec = reshape(dofConnec,size(dofConnec,1),[])';

        end

        function fInEdges = computeFinEdgesWithouRep(obj,fDisc,nodesDisc,mesh)
            s.nodesByElem = nodesDisc;
            s.type = mesh.type;
            edge = EdgesConnectivitiesComputer(s);
            edge.compute();
                        
            
            
            s.coord  = obj.coordD;
            s.connec = edge.nodesInEdges;
            s.kFace  = mesh.kFace -1;
            eM = Mesh.create(s); 


            s.edgeMesh = eM;
            s.fNodes   = fDisc;
            eF         = EdgeFunctionInterpolator(s);
            fInEdges = eF.compute()'; 

             %ehhh
             t = edge.edgesInElem;
             t2 = t; 
             t2(:,2) = t(:,3);
             t2(:,3) = t(:,2);

             t2 = t2';
             fInEdges2 = fInEdges(t2(:),:);

             fInEdges = fInEdges2;
        end


        function fEdges = computeFinEdges(obj,fDisc,nodesDisc,mesh,ndim)

            fInEdges = obj.computeFinEdgesWithouRep(fDisc,nodesDisc,mesh);
      
            %ehhh
            for iDim = 1:size(fInEdges,2)
                xc0 = fInEdges(:,iDim);
                xc = repmat(xc0',ndim*3,1);                 
                fEdges(:,iDim) = reshape(xc,[],1);
                %fEdges(:,iDim) = reshape(xc',1,[]);
            end

        end

        function fDisc = computeDiscontinousFunction(obj,fCont,connecC,connecD)
            nodesCont = reshape(connecC',1,[]);   
            nodesDisc = reshape(connecD',1,[]);   
            fDisc(nodesDisc,:) = fCont(nodesCont,:);
        end
        
        function nodesDisc = createDiscConnec(obj,connec)
            nNode  = size(connec,2);
            nElem  = size(connec,1);
            nodesDisc = 1:nNode*nElem;            
        end

        function newDofs = computeNewDofs(obj,f)
            dofsF = f.getDofConnec();
            maxDof = max(dofsF(:));
            nElem  = f.mesh.nelem;
            nEdges = f.mesh.edges.nEdgeByElem;
            nElemInEdge = 3;
            newDofsInEdge = nElemInEdge*f.ndimf;
            
            nNewDofs = nElem*nEdges*newDofsInEdge;
            newDofs  = maxDof + (1:nNewDofs);
            newDofs  = reshape(newDofs,f.ndimf,nEdges,nElemInEdge,nElem);
        end
        
    end
    
end