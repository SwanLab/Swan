classdef AssemblerFun < handle

    properties (Access = private)
        fun
        connec
    end

    methods (Access = public)
        
        function obj = AssemblerFun(cParams)
            obj.init(cParams);
        end

        function A = assemble(obj, Aelem)
            dofs   = obj.fun.computeDofConnectivity();
            nDofsE = size(dofs,1);
            nDofs  = numel(obj.fun.fValues);
            A = sparse(nDofs,nDofs);
            for i = 1:nDofsE
                for j = 1:nDofsE
                    a = squeeze(Aelem(i,j,:));
                    A = A + sparse(dofs(i,:),dofs(j,:),a,nDofs,nDofs);
                end
            end
        end

        function A = assembleFunctions(obj, Aelem, f1, f2)
            dofsF1 = f1.computeDofConnectivity();
            if isequal(f1, f2)
                dofsF2 = dofsF1;
            else
                dofsF2 = f2.computeDofConnectivity();
            end
            
            nDofs1 = numel(f1.fValues);
            nDofs2 = numel(f2.fValues);
            ndofsElem1 = size(Aelem,1);
            ndofsElem2 = size(Aelem,2);
            A = sparse(nDofs1,nDofs2);
            for i = 1:ndofsElem1
                for j = 1:ndofsElem2
                    a = squeeze(Aelem(i,j,:));
                    A = A + sparse(dofsF1(i,:),dofsF2(j,:),a,nDofs1,nDofs2);
                end
            end

        end

        function A = assembleFunctionsViaIndices(obj, Aelem, f1, f2)
            dofsF1 = f1.computeDofConnectivity()';
            if isequal(f1, f2)
                dofsF2 = dofsF1;
            else
                dofsF2 = f2.computeDofConnectivity()';
            end
            
            nElem = size(Aelem,3);
            nDofs1 = numel(f1.fValues);
            nDofs2 = numel(f2.fValues);
            ndofsElem1 = size(Aelem,1);
            ndofsElem2 = size(Aelem,2);

            res = zeros(ndofsElem1^2 * nElem, 3);
            strt = 1;
            fnsh = nElem;
            for i = 1:ndofsElem1
                dofsI = dofsF1(:,i);
                for j = 1:ndofsElem2
                    dofsJ = dofsF2(:,j);
                    a = squeeze(Aelem(i,j,:));
                    matRes = [dofsI, dofsJ, a];
                    res(strt:fnsh,:) = matRes;
                    strt = strt + nElem;
                    fnsh = fnsh + nElem;
                end
            end
            A = sparse(res(:,1), res(:,2), res(:,3), nDofs1, nDofs2);
        end

        function V = assembleV(obj, F, fun)
            % Via indices
            dofConnec = obj.fun.computeDofConnectivity();
            nDofsEl   = obj.fun.nDofsElem;
            nDofs     = obj.fun.nDofs;
            nGaus     = size(F,2);
            V = zeros(nDofs,1);
            for iDof = 1:nDofsEl
                for igaus = 1:nGaus
                    dofs = dofConnec(iDof,:);
                    c = squeeze(F(iDof,igaus,:));
                    Fadd = obj.computeAddVectorBySparse(dofs, c, nDofs);
                    % Fadd = obj.computeAddVectorByAccumarray(dofs, c, ndof);
                    V = V + Fadd;
                end
            end
        end

        function Vadd = computeAddVectorBySparse(obj,dofs, c, ndof)
           Vadd = sparse(dofs,1,c,ndof,1);
        end

        function A = assembleMatrixViaIndices(obj, Ae)
            connec    = obj.globalConnec;
%             dofConnec = obj.computeDofConnectivity()';
            nnodes  = obj.dim.nnodes;
            ndimf   = obj.dim.ndimf;
            ndofs   = ndimf*nnodes;
            nelem   = size(connec, 1);
%             nnodeEl = obj.nnodeEl;
            dofConnec = obj.computeDofConnectivity()';
            ndofEl  = size(dofConnec,2);
            res = zeros(ndofEl^2 * nelem, 3);
            strt = 1;
            fnsh = nelem;
            ndofEl1 = size(Ae,1);
            ndofEl2 = size(Ae,2);
            for i = 1:ndofEl1
                dofsI = dofConnec(:,i);
                for j = 1:ndofEl2
                    dofsJ = dofConnec(:,j);
                    a = squeeze(Ae(i,j,:));
                    matRes = [dofsI, dofsJ, a];
                    res(strt:fnsh,:) = matRes;
                    strt = strt + nelem;
                    fnsh = fnsh + nelem;
                end
            end
            A = sparse(res(:,1), res(:,2), res(:,3), ndofs, ndofs);
        end

        function F = assembleVectorStokes(obj, FelemCell, f1, f2)
            % Stokes
            funs = {f1,f2};
            nFields = numel(funs);
            b_global = cell(nFields,1);
            for iField = 1:nFields
                f = funs{iField};
                Felem = FelemCell{iField,1};
                dofsElem = f.computeDofConnectivity();
                nDofs  = numel(f.fValues);
                nDofsE = size(Felem,1);
                nGaus  = size(Felem,2);
                b = zeros(nDofs,1);
                for iDof = 1:nDofsE
                    for iGaus = 1:nGaus
                        c = squeeze(Felem(iDof,iGaus,:));
                        idof_elem = dofsElem(iDof,:);
                        b = b + sparse(idof_elem,1,c',nDofs,1);
                    end
                end
                b_global{iField,1} = b;
            end
            F = cell2mat(b_global);

        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.fun    = cParams.fun;
        end

    end

end