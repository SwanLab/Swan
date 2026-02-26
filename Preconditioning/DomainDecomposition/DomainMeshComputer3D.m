classdef DomainMeshComputer3D < handle

    properties (Access = public)
        domainMesh
        domainMeshDisc
        localGlobalConnec
    end

    properties (Access = private)
        meshReference
        nSubdomains
        interfaceMeshSubDomain
        ninterfaces
        meshSubDomain
        interfaceConnec
        tolSameNode
    end

    properties (Access = private)
        connecGlob
        coordGlob
        updtConnecGlob
        updtCoordGlob
    end

    methods (Access = public)

        function obj = DomainMeshComputer3D(cParams)
            obj.init(cParams)

        end

        function compute(obj)
            obj.connecGlob = obj.createGlobalConnec(obj.meshReference);
            obj.createGlobalCoord();
            obj.updateGlobalConnec();
            obj.updateGlobalCoord();
            obj.computeLocalGlobalConnec()
            s.coord        = obj.updtCoordGlob;
            s.connec       = obj.updtConnecGlob;
            obj.domainMesh = Mesh.create(s);

            s.connec = obj.connecGlob;
            s.coord = obj.coordGlob;
            obj.domainMeshDisc = Mesh.create(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.tolSameNode     = cParams.tolSameNode;
            obj.meshReference   = cParams.meshReference;
            obj.nSubdomains     = cParams.nSubdomains;
            obj.interfaceMeshSubDomain = cParams.interfaceMeshSubDomain;
            obj.ninterfaces     = cParams.ninterfaces;
            obj.meshSubDomain   = cParams.meshSubDomain;
            obj.interfaceConnec = cParams.interfaceConnec;
        end

        function gConnec = createGlobalConnec(obj,m)
            % we simply create a connectivity matrix as if no nodes are
            % shared, just summing nnodes. We need that to create the
            % global connectivity matrix that considers shared nodes.
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            nZ = obj.nSubdomains(3);
            nelem      = m.nelem;
            nnodeElem  = m.nnodeElem;
            nnodes     = m.nnodes;
            connec0    = m.connec;
            %             dConnec = connec0 + nnodes*(nX*(jDom-1)+iDom-1);
            gConnec=zeros(nX*nY*nZ*nelem,nnodeElem);
            for kDom = 1:nZ
                for jDom = 1:nY
                    for iDom = 1:nX
                        indLinear = nX*nY*(kDom-1)+nX*(jDom-1)+iDom;
                        rowIn     = (indLinear-1)*nelem+1;
                        rowEnd    = indLinear*nelem;
                        gConnec(rowIn:rowEnd,:)=connec0+nnodes*(indLinear-1);
                    end
                end
            end
        end

        function  createGlobalCoord(obj)
            % we simply create the coordinate matrix as if no nodes are
            % shared, just summing nnodes.
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            nZ = obj.nSubdomains(3);
            ndim     = obj.meshReference.ndim;
            nnodes  = obj.meshReference.nnodes;
            meshsd = obj.meshSubDomain;
            %             dConnec = connec0 + nnodes*(nX*(jDom-1)+iDom-1);
            coordGlob=zeros(nX*nY*nZ*nnodes,ndim);
            for kDom = 1:nZ
                for jDom = 1:nY
                    for iDom = 1:nX
                        indLinear= nX*nY*(kDom-1)+nX*(jDom-1)+iDom;
                        rowIn=(indLinear-1)*nnodes+1;
                        rowEnd=indLinear*nnodes;
                        coordGlob(rowIn:rowEnd,:)=meshsd{jDom,iDom,kDom}.coord;
                    end
                end
            end
            obj.coordGlob=coordGlob;
        end

        function   updateGlobalConnec(obj)
            % I leave this here to know what i was doing. Below, the
            % optimized version from chatgpt
% %             updtConnecGlob = obj.connecGlob;
% %             rCconnec       = obj.interfaceConnec;
% %             nref           = size(rCconnec,1);
% %             %             ncopies        = size(rCconnec,2);
% %             for iref=1:nref
% %                 aux = rCconnec(iref,:);
% %                 aux = aux(aux>0);
% %                 for icopy=2:length(aux)
% %                     updtConnecGlob(updtConnecGlob==aux(icopy))  = aux(1);
% %                     updtConnecGlob(updtConnecGlob>aux(icopy))   = updtConnecGlob(updtConnecGlob>aux(icopy))-1;
% %                     aux(aux>aux(icopy)) =  aux(aux>aux(icopy))-1;
% %                     rCconnec(rCconnec>aux(icopy)) = rCconnec(rCconnec>aux(icopy))-1;
% %                 end
% %             end
% %             obj.updtConnecGlob=updtConnecGlob;

            updtConnecGlob = obj.connecGlob;
            rCconnec       = obj.interfaceConnec;

            % Step 1: Build map from duplicate → master
            duplicateMap = containers.Map('KeyType', 'uint32', 'ValueType', 'uint32');

            for i = 1:size(rCconnec, 1)
                group = rCconnec(i, :);
                group = group(group > 0);
                master = group(1);
                for j = 2:length(group)
                    duplicateMap(group(j)) = master;
                end
            end

            % Step 2: Replace duplicates in connectivity
            allKeys = cell2mat(duplicateMap.keys);
            allVals = cell2mat(duplicateMap.values);

            % Flatten connectivity to vector
            originalSize = size(updtConnecGlob);
            flatConnec = updtConnecGlob(:);

            % Replace duplicates
            [found, idx] = ismember(flatConnec, allKeys);
            flatConnec(found) = allVals(idx(found));

            % Step 3: Compact node indices starting from 1
            [~, ~, compactConnec] = unique(flatConnec);

            % Reshape back to original connectivity matrix size
            updtConnecGlob = reshape(compactConnec, originalSize);

            % Output:
            obj.updtConnecGlob = updtConnecGlob;
        end

        function  updateGlobalCoord(obj)
            tol = obj.tolSameNode;
            nodes = obj.connecGlob(:);
            nodes = unique(nodes);
            for i = 1:size(obj.interfaceConnec,2)
                [~,ind] = ismember(obj.interfaceConnec(:,i),nodes);
                nodes(ind(ind>0)) = obj.interfaceConnec(ind>0,1);
            end
            %             [~,ind] = ismember(obj.interfaceConnec(:,2),nodes);
            %             nodes(ind) = obj.interfaceConnec(:,1);
            nodes= unique(nodes,'stable');
            uniqueVals = obj.coordGlob(nodes,:);
            obj.updtCoordGlob = uniqueVals;
            %             A = obj.coordGlob;
            %
            %             % Initialize a logical index to track unique rows
            %             %tic
            %             isUnique = true(size(A, 1), 1);
            %
            %             for i = 1:size(A, 1)
            %                 if isUnique(i)
            %                     % Compute row-wise tolerance checks
            %                     rowDiff = vecnorm(A - A(i, :),2,2);
            %                     isSimilarRow = all(rowDiff <= tol, 2); % Rows similar within tolerance
            %                     isUnique = isUnique & ~isSimilarRow;  % Mark similar rows as non-unique
            %                     isUnique(i) = true;                  % Keep the first occurrence
            %                 end
            %             end
            %
            %             uniqueVals = A(isUnique,:);
            %             % toc
            %             %
            %             % tic
            %             % % Compute pairwise differences between rows
            %             % diffMatrix = abs(A - permute(A, [3, 2, 1]));  % Compute differences between rows
            %             % maxDiff = squeeze(max(diffMatrix, [], 2));    % Maximum difference across columns
            %             %
            %             % % Find unique rows by checking if they are within the tolerance
            %             % isSimilar = maxDiff <= tol;                   % Rows within tolerance
            %             % isUnique = ~any(tril(isSimilar, -1), 2);      % Ignore duplicates in lower triangle
            %             %
            %             % % Extract unique rows
            %             % uniqueRows = A(isUnique, :);
            %             % toc
            %
            %             obj.updtCoordGlob = uniqueVals;
            %             %tol = obj.tolSameNode;
            %           %  [~, colindices] = uniquetol(obj.coordGlob,tol, 'ByRows', true);   %get indices of unique value. Is sorted BY VALUE
            %          %   obj.updtCoordGlob = obj.coordGlob(sort(colindices),:);
            %
            %             %obj.updtCoordGlob  = unique(obj.coordGlob,'rows','stable');%'stable');
        end

        function  computeLocalGlobalConnec(obj)
            nodeG  = reshape(obj.updtConnecGlob',[],1);
            nodeL  = reshape(obj.connecGlob',[],1);
            nnode = obj.meshReference.nnodes;
            localGlobalConnec(:,1) = nodeG;
            localGlobalConnec(:,2) = nodeL;
            localGlobalConnec      = unique(localGlobalConnec,'rows','stable');
            localGlobalConnec      = permute(reshape(localGlobalConnec',2,nnode,[]),[2,1,3]);

            for dom = 1:size(localGlobalConnec,3)
                localGlobalConnec(:,2,dom) = localGlobalConnec(:,2,dom)-(dom-1)*nnode;
            end
            obj.localGlobalConnec  = localGlobalConnec;
        end
    end
end