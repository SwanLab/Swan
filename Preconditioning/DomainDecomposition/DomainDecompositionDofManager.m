classdef DomainDecompositionDofManager < handle

    properties (Access = public)
        interfaceConnec
        interfaceConnecReshaped
        interfaceDom
        intConecLocal
    end

    properties (Access = private)
        nDof
        nReferenceDof
        localGlobalDof
        interfaceDof
        CachedLocalRows   
        CachedColIdx
        CachedGlobalDest
    end

    properties (Access = private)
        nSubdomains
        locGlobConnec
        nBoundaryNodes
        nReferenceNodes
        nDimf
        nNodes
    end

    methods (Access = public)

        function init(obj,cParams)
            obj.nSubdomains     = cParams.nSubdomains;
            obj.interfaceConnec = cParams.interfaceConnec;
            obj.interfaceConnecReshaped = cParams.interfaceConnecReshaped;
            obj.locGlobConnec   = cParams.locGlobConnec;
            obj.nBoundaryNodes  = cParams.nBoundaryNodes;
            obj.nReferenceNodes = cParams.nReferenceNodes;
            obj.nNodes          = cParams.nNodes;
            obj.nDimf           = cParams.nDimf;
            obj.nDof            = obj.nNodes*obj.nDimf;
            obj.nReferenceDof   = obj.nReferenceNodes*obj.nDimf;
        end

        function obj = DomainDecompositionDofManager(cParams)
            obj.init(cParams)
            obj.createlocalGlobalDofConnec();
            obj.computeLocalInterfaceDof();
            obj.SetupIndices();
        end

        %         function f = scaleInterfaceValues(obj,f,w)
        %             nint = numel(obj.interfaceDof);
        %             weight = [w,1-w];
        %             for iint = 1:nint
        %                 dofI = obj.interfaceDof{iint};
        %                 ndom = size(dofI,2);
        %                 for idom = 1:ndom
        %                     dom = obj.interfaceDom(iint,idom);
        %                     dof = dofI(:,idom);
        %                     f(dof,dom) = weight(idom)* f(dof,dom);
        %                 end
        %             end
        %         end

%         function f = scaleInterfaceValues(obj, f, w)
%             % Flatten data structures to avoid nested loops
%             allDofs = cell2mat(obj.interfaceDof(:));
%             nRowsPerCell = cellfun(@(x) size(x, 1), obj.interfaceDof(:));
%             allDoms = repelem(obj.interfaceDom, nRowsPerCell, 1);
% 
%             weights = [w, 1-w];
%             nrows_f = size(f, 1);
% 
%             % Process each column (domain side) vectorially
%             for idom = 1:size(allDofs, 2)
%                 rows = allDofs(:, idom);
%                 cols = allDoms(:, idom);
% 
%                 linIdx = rows + (cols - 1) * nrows_f;
%                 f(linIdx) = f(linIdx) * weights(idom);
%             end
%         end

        function f = scaleInterfaceValues(obj, f, w)

    % Flatten DOFs
    allDofs = cell2mat(obj.interfaceDof(:));   % (N × 2)

    % Compute number of rows per cell using a fast loop
    c = obj.interfaceDof(:);
    n = numel(c);
    nRowsPerCell = zeros(n,1);
    for i = 1:n
        nRowsPerCell(i) = size(c{i},1);
    end

    % Expand domain indices
    allDoms = repelem(obj.interfaceDom, nRowsPerCell, 1);  % (N × 2)

    % Precompute weights
    weights = [w, 1-w];

    % Compute linear indices for all entries at once
    nrows_f = size(f,1);
    linIdx = allDofs + (allDoms - 1) * nrows_f;  % (N × 2)

    % Apply scaling vectorially
    f(linIdx) = f(linIdx) .* weights;

        end

%         function f = scaleInterfaceValues(obj, f, w)
%     % 1. Get Cell Content (Vectorized)
%     % Flatten to column to ensure consistent ordering
%     dofCells = obj.interfaceDof(:);
%     
%     % 2. Count Rows (Replaces the loop)
%     % cellfun with 'size' is extremely fast and internal
%     nRowsPerCell = cellfun('size', dofCells, 1);
% 
%     % 3. Flatten DOFs (Replaces cell2mat)
%     % vertcat is much faster than cell2mat for matrices
%     allDofs = vertcat(dofCells{:});
% 
%     % 4. Expand Domains
%     % Assumes interfaceDom matches the order of interfaceDof
%     allDoms = repelem(obj.interfaceDom, nRowsPerCell, 1);
% 
%     % 5. Linear Indices & Scaling
%     stride = size(f, 1);
%     linIdx = allDofs + (allDoms - 1) * stride;
%     
%     % 6. Apply Weights
%     % Implicit expansion handles the [w, 1-w] broadcasting automatically
%     f(linIdx) = f(linIdx) .* [w, 1-w];
% end


        function m = scaleInterfaceValuesMatrix(obj,m,w)
            nint = numel(obj.interfaceDof);
            weight = [w,1-w];
            for iint = 1:nint
                dofI = obj.interfaceDof{iint};
                ndom = size(dofI,2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    dof = dofI(:,idom);
                    m(dof,dof,dom) = weight(idom)* m(dof,dof,dom);
                end
            end
        end

        function fG = local2global(obj,fL)
            fG   = zeros(obj.nDof,obj.nSubdomains(1)*obj.nSubdomains(2));
            ind    = 1;
            for jdom = 1: obj.nSubdomains(2)
                for idom = 1: obj.nSubdomains(1)
                    lDof  = obj.localGlobalDof{jdom,idom}(:,1);
                    gDof = obj.localGlobalDof{jdom,idom}(:,2);
                    fG(lDof,ind) = fL(gDof,ind);
                    ind=ind+1;
                end
            end
        end
        %
        %         function fG = AssembleLocal2GlobalVector(obj,fL)
        %             fG   = zeros(obj.nDof,1);
        %             ind    = 1;
        %             for jdom = 1: obj.nSubdomains(2)
        %                 for idom = 1: obj.nSubdomains(1)
        %                     lDof  = obj.localGlobalDof{jdom,idom}(:,1);
        %                     gDof = obj.localGlobalDof{jdom,idom}(:,2);
        %                     fG(lDof) = fG(lDof)+ fL(gDof,ind);
        %                     ind=ind+1;
        %                 end
        %             end
        %         end

%         function fG = AssembleLocal2GlobalVector(obj, fL)
%             mapCell = obj.localGlobalDof.';
%             bigMap = cell2mat(mapCell(:));
% 
%             global_dest_idx = bigMap(:, 1);
%             local_row_idx   = bigMap(:, 2);
% 
%             nEntriesPerSubdomain = size(mapCell{1}, 1);
%             nSubdomains          = numel(mapCell);
% 
%             col_idx = repelem(1:nSubdomains, nEntriesPerSubdomain).';
% 
%             idx_fL = local_row_idx + (col_idx - 1) * size(fL, 1);
%             values = fL(idx_fL);
% 
%             fG = accumarray(global_dest_idx, values, [obj.nDof, 1]);
%         end

function fG = AssembleLocal2GlobalVector(obj, fL)
    stride = size(fL, 1);
    read_idx = obj.CachedLocalRows + (obj.CachedColIdx * uint32(stride));
    values = fL(read_idx);
    fG = accumarray(obj.CachedGlobalDest, values, [obj.nDof, 1]);
end




%         function fL = global2local(obj,fG)
%             ndimf  = obj.nDimf;
%             fL     = zeros(obj.nReferenceNodes*ndimf,obj.nSubdomains(1)*obj.nSubdomains(2));
%             ind    = 1;
%             for jdom = 1:obj.nSubdomains(2)
%                 for idom = 1:obj.nSubdomains(1)
%                     lGconnec = obj.localGlobalDof{jdom,idom};
%                     fL(lGconnec(:,2),ind) = fG(lGconnec(:,1));
%                     ind=ind+1;
%                 end
%             end
%         end

        function fL = global2local(obj, fG)
    % 1. Initialize fL
    % We calculate the stride (rows per element) dynamically
    stride = obj.nReferenceNodes * obj.nDimf;
    nTotalElements = prod(obj.nSubdomains);
    
    fL = zeros(stride, nTotalElements);

    % 2. Calculate Linear Write Indices (Where to put data in fL)
    % We reuse the cached properties from SetupIndices
    % Index = Local_Row + (Column_Index * Stride)
    write_idx = obj.CachedLocalRows + (obj.CachedColIdx * uint32(stride));

    % 3. Scatter Data (Zero Loops)
    % We read from Global (fG) and write to Local (fL)
    % MATLAB automatically handles the "fan-out" (copying one global value 
    % to multiple local elements that share that node).
    fL(write_idx) = fG(obj.CachedGlobalDest);
end

        function mL = global2localMatrix(obj,mG)
            ndimf  = obj.nDimf;
            mL     = zeros(obj.nReferenceNodes*ndimf,obj.nReferenceNodes*ndimf,obj.nSubdomains(1)*obj.nSubdomains(2));
            ind    = 1;
            for jdom = 1:obj.nSubdomains(2)
                for idom = 1:obj.nSubdomains(1)
                    lGconnec = obj.localGlobalDof{jdom,idom};
                    mL(lGconnec(:,2),lGconnec(:,2),ind) = mG(lGconnec(:,1),lGconnec(:,1));
                    ind=ind+1;
                end
            end
        end

    end

    methods (Access = private)

   function SetupIndices(obj)
    % 1. Enforce Row-Major Order (Critical!)
    % Your original code transposed the cell array. We must do the same 
    % so that Subdomain 1 corresponds to fL column 1, etc.
    orderedCells = obj.localGlobalDof.'; 
    
    % 2. Get sizes of all 30 subdomains (Robust to variable sizes)
    % orderedCells is now 10x3 (linearized as 30x1)
    counts = cellfun('size', orderedCells, 1);
    
    % 3. Flatten the Map
    % efficient stacking of all 30 matrices into one big list
    bigMap = vertcat(orderedCells{:});

    % 4. Generate Column Indices
    % If Subdomain 1 has 896 entries, we repeat "Col 0" 896 times.
    % If Subdomain 2 has 896 entries, we repeat "Col 1" 896 times.
    nSubdomains = numel(orderedCells);
    col_indices = repelem(0:nSubdomains-1, counts(:)); 

    % --- STORE PROPERTIES (uint32 for speed/memory) ---
    
    % Local Row Indices (Column 2 of your map)
    obj.CachedLocalRows = uint32(bigMap(:, 2));
    
    % Subdomain Column Indices (0, 0... 1, 1... 29, 29...)
    obj.CachedColIdx    = uint32(col_indices(:)); 
    
    % Global Destination Indices (Column 1 of your map)
    obj.CachedGlobalDest = uint32(bigMap(:, 1));
end
        function  createlocalGlobalDofConnec(obj)
            ndimf = obj.nDimf;
            ndom  = obj.nSubdomains(1)*obj.nSubdomains(2);
            for dom = 1:ndom
                row = ceil(dom/obj.nSubdomains(1));
                col = dom-(row-1)*obj.nSubdomains(1);
                localGlobalDofConnecDom = zeros(1,2);
                nodeG = obj.locGlobConnec(:,1,dom);
                nodeL = obj.locGlobConnec(:,2,dom);
                for iunkn = 1:ndimf
                    dofConec = [ndimf*(nodeG - 1) + iunkn ,  ndimf*(nodeL - 1) + iunkn] ;
                    localGlobalDofConnecDom = [localGlobalDofConnecDom;dofConec];
                end
                localGlobalDofConnec{row,col} = localGlobalDofConnecDom(2:end,:);
            end
            obj.localGlobalDof = localGlobalDofConnec;
        end

        function computeLocalInterfaceDof(obj)

            intConecResh = obj.interfaceConnecReshaped;



            nint = numel(intConecResh);
            ndimf = obj.nDimf;
            ndofs = obj.nReferenceNodes*ndimf;
            for iint=1:nint
                intConecIint = intConecResh{iint};
                ndom = size(intConecIint,2); %length(intConec(1,:,iint));
                intConecL = zeros(size(intConecIint));
                interfaceDof = [];
                for idom = 1:ndom
                    dofaux=0;

                    nodesI = intConecIint(:,idom);
                    dom = ceil(intConecIint(1,idom)/obj.nReferenceNodes);
                    globaldof = (dom-1)*ndofs;
                    for iunkn=1:ndimf
                        DOF = ndimf*(nodesI - 1) + iunkn;
                        DOF = DOF-globaldof;
                        dofaux= [dofaux; DOF];
                    end
                    interfaceDof(:,idom) = dofaux(2:end);
                    interfaceDom(iint,idom) = dom;
                    intConecL(:,idom) = nodesI - (dom-1)*obj.nReferenceNodes;
                end
                intConecLiint{iint} = intConecL(:,idom);
                interfaceDofIint{iint} = interfaceDof;
            end
            obj.interfaceDof  = interfaceDofIint;
            obj.interfaceDom  = interfaceDom;
            obj.intConecLocal = intConecLiint;

        end

    end

end