classdef P1Refiner < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = P1Refiner(fP1)
            obj.computeDofConnec(fP1)
            
        end
        
    end
    
    methods (Access = private)
        
        function computeDofConnec(obj,f)
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

            connec(:,1,:) = [nV1 ; e1d1 ;e3d3];
            connec(:,2,:) = [e1d3; nV2 ;e2d1];
            connec(:,3,:) = [e1d2; e2d2; e3d2];
            connec(:,4,:) = [e2d3;  nV3; e3d1];
            %connec = permute(connec,[2,1,3]);
            newConnec = reshape(connec,size(connec,1),[])';

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