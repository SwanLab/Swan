classdef TetradecahedronMeshComputer < handle

    properties (Access = private, Constant)
        planeCoeffs = [1  0  0   1     1    1     1;
            0  1  0   1     1   -1    -1;
            0  0  1   1    -1    1    -1;
            0  0  0  9/4   5/4  5/4   1/4;
            1  1  1  3/4  -1/4 -1/4  -5/4];
        center = [0.5,0.5,0.5];
    end

    properties (Access = private)
        mesh
        boundaryNodes
    end

    methods (Access = public)

        function obj = TetradecahedronMeshComputer(fileName)
            a.fileName = fileName;
            fem = FemDataContainer(a);
            obj.mesh = fem.mesh;
            obj.boundaryNodes = fem.newBC.microSlave(:,1);
        end

        function mesh = getMesh(obj)
            mesh = obj.mesh;
        end

        function MS = getMasterSlave(obj)
            MSfull = obj.computeMasterSlaveFull();
            MS     = obj.reduceMasterSlave(MSfull);
        end

    end

    methods (Access = private)

        function MSfull = computeMasterSlaveFull(obj)
            nodes = 1:obj.mesh.nnodes;
            coordNodes = [nodes',obj.mesh.coord];
            coordFaces = coordNodes(obj.boundaryNodes,:);

            MSfull = [];
            for i=1:size(obj.planeCoeffs,2)
                n = obj.planeCoeffs(1:3,i);

                fMaster = (abs(coordFaces(:,2:4)*n - obj.planeCoeffs(4,i)) < 1e-2);
                fSlave  = (abs(coordFaces(:,2:4)*n - obj.planeCoeffs(5,i)) < 1e-2);
                coordMaster = coordFaces(fMaster,:);
                coordSlave  = coordFaces(fSlave,:);

                distVec = obj.computeDistanceVector(coordMaster,n);
                coordMasterNew = [coordMaster(:,1), coordMaster(:,2:end) + distVec];
                cM = sortrows(coordMasterNew,[2:4]); cS = sortrows(coordSlave,[2:4]);
                %max(abs(cM(:,2:4)-cS(:,2:4)),[],'all') % Check position difference between nodes
                MSfull = [MSfull; [cM(:,1) cS(:,1)]];
            end
        end

        function [distVec] = computeDistanceVector(obj,coordFace,n)
            centerVec = coordFace(1,2:4) - obj.center;
            distVec = -2*(dot(n,centerVec)/norm(n)^2)*n';
        end

        function MSreduced = reduceMasterSlave(obj,MSfull)
            % Remove edges and vertices
            vals = MSfull(:);
            [uv,~,idx] = unique(vals);
            counts = accumarray(idx,1);
            repeatedVals = uv(counts == 1);
            rowsToKeep = any(ismember(MSfull,repeatedVals),2);
            MSfaces = MSfull(rowsToKeep,:);
            % Obtain edges
            edges = uv(counts == 2);
            rowsToKeep = any(ismember(MSfull,edges),2);
            MSedges = MSfull(rowsToKeep,:);
            MSedges = obj.removeAndOrder(MSedges);
            %Obtain vertices
            edges = uv(counts == 3);
            rowsToKeep = any(ismember(MSfull,edges),2);
            MSvertices = MSfull(rowsToKeep,:);
            MSvertices = obj.removeAndOrder(MSvertices);

            MSreduced = [MSfaces;MSedges;MSvertices];
        end

        function MSnew = removeAndOrder(~,MS)
            MSnew = [];
            while(~isempty(MS))
                idx2add = find(MS(:,1) == MS(1,1));
                idx2permute = find(MS(:,2) == MS(1,1));
                MSnew = [MSnew; MS(idx2add,:); [MS(idx2permute,2), MS(idx2permute,1)]];

                unique(MS([idx2add;idx2permute],:));
                vals2remove = unique(MS([idx2add;idx2permute],:));

                idx2remove = any(ismember(MS,vals2remove),2);
                MS(idx2remove,:) = [];
            end
        end

    end

end
