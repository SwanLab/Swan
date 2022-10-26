classdef Assembler < handle

    properties (Access = private)
        dim
        globalConnec
        nnodeEl
    end

    methods (Access = public)
        
        function obj = Assembler(cParams)
            obj.init(cParams);
        end

        function A = assemble(obj, aElem)
%             disp('normal')
%             tic
%                 A = obj.assembleMatrix(aElem);
%             toc
%             disp('nou')
%             tic
                A = obj.assembleMatrixViaIndices(aElem);
%             toc
        end

        function B = assembleB(obj, Bfull, d)
            B = obj.assembleBMatrix(Bfull, d);
        end

        function C = assembleC(obj, Cmat, dvol)
            C = obj.assembleCMatrix(Cmat, dvol);
        end

        function V = assembleV(obj, F)
            V = obj.assembleVector(F);
        end

        function A = assembleFields(obj, Ae, f1, f2)
            switch f1.galerkinType
                case 'CONTINUOUS'
                    A = obj.assembleFieldsCG(Ae, f1, f2);
                case 'DISCONTINUOUS'
                    A = obj.assembleFieldsDG(Ae, f1, f2);
            end
        end

        function A = assembleVectorFields(obj, Ae, f1, f2)
            A = obj.assembleVectorWithFields(Ae, f1, f2);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.dim          = cParams.dim;
            obj.nnodeEl      = cParams.nnodeEl;
            obj.globalConnec = cParams.globalConnec;
        end

        function A = assembleMatrix(obj, Ae)
            connec    = obj.globalConnec;
            dofConnec = obj.computeDofConnectivity()';
            ndofs   = obj.dim.ndofs;
            ndimf   = obj.dim.ndimf;
            nnodeEl = obj.nnodeEl;
            A = sparse(ndofs,ndofs);
            for i = 1:nnodeEl*ndimf
                dofsI = dofConnec(:,i);
                for j = 1:nnodeEl*ndimf
                    dofsJ = dofConnec(:,j);
                    a = squeeze(Ae(i,j,:));
                    Aadd = obj.computeAaddBySparse(a, dofsI, dofsJ);
%                     Aadd = obj.computeAaddByAccumarray(a, dofsI, dofsJ);
                    A = A + Aadd;
                end
            end
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
        
        function dofConnec = computeDofConnectivity(obj)
            connec  = obj.globalConnec;
            ndimf   = obj.dim.ndimf;
            nnodeEl = size(connec, 2); % obj.dim.nnodeElem
            ndofsEl = nnodeEl * ndimf; %obj.dim.ndofsElem;
            dofsElem  = zeros(ndofsEl,size(connec,1));
            for inode = 1:nnodeEl
                for iunkn = 1:ndimf
                    idofElem   = obj.nod2dof(inode,iunkn);
                    globalNode = connec(:,inode);
                    idofGlobal = obj.nod2dof(globalNode,iunkn);
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            dofConnec = dofsElem;
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimf;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

        function Cadd = computeAaddBySparse(obj,a, dofsI, dofsJ)
           d = obj.dim;
           ndofs = d.ndofs;
           Cadd = sparse(dofsI,dofsJ,a,ndofs,ndofs);
        end

        function Cadd = computeAaddByAccumarray(obj,a, dofsI, dofsJ)
           d = obj.dim;
           ndofs = d.ndofs;
           index = [dofsI, dofsJ];
           Cadd = accumarray(index,a,[ndofs ndofs],[],[],true);
        end

        %% Vector assembly

        function V = assembleVector(obj, F)
            dofsInElem = obj.computeDofConnectivity();
            ndofPerElem = obj.dim.ndofsElem;
            ndof        = obj.dim.ndofs;
            V = zeros(ndof,1);
            for iDof = 1:ndofPerElem
                dofs = dofsInElem(iDof,:);
                c = F(iDof,:);
                Fadd = obj.computeAddVectorBySparse(dofs, c, ndof);
                % Fadd = obj.computeAddVectorByAccumarray(dofs, c, ndof);
                V = V + Fadd;
            end
        end

        function Vadd = computeAddVectorBySparse(obj,dofs, c, ndof)
           Vadd = sparse(dofs,1,c',ndof,1);
        end

        function Vadd = computeAddVectorByAccumarray(obj,dofs,c, ndof)
           Vadd = accumarray(dofs',c',[ndof 1]);
        end
        
        %% B matrix assembly
        function Bt = assembleBMatrix(obj, Bfull, dims)
            d = obj.dim;
            ntot  = size(Bfull,1);
            dims.nvoigt = ntot/(dims.nelem * dims.ngaus);
            ndofGlob = d.ndof;
            Bt = sparse(ntot,ndofGlob);
            for idof = 1:d.ndofPerElement
                dofs  = obj.computeGlobalDofs(idof, dims);
                Bidof = Bfull(:,idof);
                Bdof = obj.computeBdofBySparse(Bidof,dofs);
%                 Bdof = obj.computeBdofByAccumarray(Bidof,dofs);
                Bt = Bt + Bdof;
            end
        end

        function dofs = computeGlobalDofs(obj, idof, dims)
            nvoigt = dims.nvoigt;
            gDofs = obj.transformLocal2Global(idof);
            dofs = repmat(gDofs',dims.ngaus*nvoigt,1);
            dofs = dofs(:);
        end

        function gDofs = transformLocal2Global(obj,iDof)
            d       = obj.dim;
            ndimf   = d.ndimf;
            nnodeEl = obj.nnodeEl;
            nodes        = obj.globalConnec;
            nodesInElem  = reshape(repmat(1:nnodeEl,ndimf,1),1,[]);
            dofs         = repmat(1:ndimf,1,nnodeEl);
            inode        = nodesInElem(iDof);
            iunkn        = dofs(iDof);
            nodeI        = nodes(:,inode);
            gDofs   = ndimf*(nodeI-1) + iunkn;
        end

        function Bdof = computeBdofBySparse(obj,Bidof,dofs)
            d = obj.dim;
            ntot  = size(Bidof,1);
            ndofGlob = d.ndof;
            Bdof = sparse(1:ntot,dofs,Bidof,ntot,ndofGlob);
        end

        function Bdof = computeBdofByAccumarray(obj,Bidof,dofs)
            d = obj.dim;
            ntot  = size(Bidof,1);
            ndof = d.ndof;
            posI = 1:ntot;
            index = [posI', dofs];
%             Bdof = accumarray(dofs,Bidof,[ntot 1]);
            Bdof = accumarray(index,Bidof,[ntot ndof],[],[],true);
        end

        %% C matrix assembly

       function CmatTot = assembleCMatrix(obj, Cmat, dvol)
           nvoigt = size(Cmat,1);
           nelem  = size(Cmat,3);
           ngaus  = size(dvol,2);
           ntot   = ngaus*nelem*nvoigt;
           CmatTot = sparse(ntot,ntot);
           for istre = 1:nvoigt
               for jstre = 1:nvoigt
                   for igaus = 1:ngaus
                       posI = (istre)+(nvoigt)*(igaus-1) : ngaus*nvoigt : ntot;
                       posJ = (jstre)+(nvoigt)*(igaus-1) : ngaus*nvoigt : ntot;
                       Ci = Cmat(istre,jstre,:,igaus);
                       Ct = squeeze(Ci).*dvol(:,igaus);
                       Cadd = obj.computeCaddBySparse(Ct, posI, posJ, ntot);
                       CmatTot = CmatTot + Cadd;
                   end
               end
           end
       end

       function Cadd = computeCaddBySparse(obj,Ct, posI, posJ, ntot)
           Cadd = sparse(posI,posJ,Ct,ntot,ntot);
       end

       function Cadd = computeCaddByAccumarray(obj,Ct, posI, posJ, ntot)
           index = [posI', posJ'];
           Cadd = accumarray(index,Ct,[ntot ntot],[],[],true);
       end


        %% With Fields
      
        function A = assembleFieldsCG(obj, Aelem, f1, f2)
            % CONTINUOUS GALERKIN
            % Can be accelerated using indices
            dofsF1 = obj.computeFieldDofs(f1);
            if isequal(f1, f2)
                dofsF2 = dofsF1;
            else
                dofsF2 = obj.computeFieldDofs(f2);
            end
            
            ndofs1 = f1.dim.ndofs;
            ndofs2 = f2.dim.ndofs;
            ndofsElem1 = f1.dim.ndofsElem;
            ndofsElem2 = f2.dim.ndofsElem;
            A = sparse(ndofs1,ndofs2);
            for i = 1:ndofsElem1
                for j = 1:ndofsElem2
                    a = squeeze(Aelem(i,j,:));
                    A = A + sparse(dofsF1(i,:),dofsF2(j,:),a,ndofs1,ndofs2);
                end
            end
        end
      
        function A = assembleFieldsDG(obj, Aelem, f1, f2)
            % DISCONTINUOUS GALERKIN
            % Can be accelerated using indices
            dofs = 1:f1.dim.ndofs;
            nElem = size(Aelem,3);
            dofsF1 = reshape(dofs, [f1.dim.nnodeElem,nElem]);
            dofsF2 = dofsF1;
            
            ndofs1 = f1.dim.ndofs;
            ndofs2 = f2.dim.ndofs;
            ndofsElem1 = f1.dim.ndofsElem;
            ndofsElem2 = f2.dim.ndofsElem;
            A = sparse(ndofs1,ndofs2);
            for i = 1:ndofsElem1
                for j = 1:ndofsElem2
                    a = squeeze(Aelem(i,j,:));
                    A = A + sparse(dofsF1(i,:),dofsF2(j,:),a,ndofs1,ndofs2);
                end
            end
        end
        
      
        function F = assembleVectorWithFields(obj, FelemCell, f1, f2)
            fields = {f1,f2};
            nfields = numel(fields);
            for ifield = 1:nfields
                field = fields{ifield};
                dims  = field.dim;
                Felem = FelemCell{ifield,1};
                dofsElem = obj.computeFieldDofs(field);
                b = zeros(dims.ndofs,1);
                for i = 1:dims.nnodeElem*dims.ndimf
                    for igaus = 1:size(Felem,2)
                    c = squeeze(Felem(i,igaus,:));
                    idof_elem = dofsElem(i,:);
                    b = b + sparse(idof_elem,1,c',dims.ndofs,1);
                    end
                end
                b_global{ifield,1} = b;
            end
            F =cell2mat(b_global);

        end

        function dofConnec = computeFieldDofs(obj, field)
            connec  = field.connec;
            ndimf   = field.dim.ndimf;
            nnodeEl = size(connec, 2); % obj.dim.nnodeElem
            ndofsEl = nnodeEl * ndimf; %obj.dim.ndofsElem;
            dofsElem  = zeros(ndofsEl,size(connec,1));
%             idof(:,1)= ndimf*(inode - 1) + iunkn;
            for inode = 1:nnodeEl
                for iunkn = 1:ndimf
                    idofElem   = ndimf*(inode - 1) + iunkn;
                    globalNode = connec(:,inode);
                    idofGlobal = ndimf*(globalNode - 1) + iunkn;
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            dofConnec = dofsElem;
        end

    end

end