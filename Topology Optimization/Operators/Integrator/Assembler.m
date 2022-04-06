classdef Assembler < handle

    properties (Access = private)
        dim
        globalConnec
        %dofsInElem
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

        function B = assembleB(obj, Bfull)
            B = obj.assembleBMatrix(Bfull);
        end

        function C = assembleC(obj, Cmat, dvol)
            C = obj.assembleCMatrix(Cmat, dvol);
        end

        function V = assembleV(obj, F, dofs)
            V = obj.assembleVector(F, dofs);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.dim          = cParams.dim;
            obj.globalConnec = cParams.globalConnec;
        end

        function A = assembleMatrix(obj, aElem)
            dofConnec = obj.computeDofConnectivity()';
            ndofs  = obj.dim.ndof;
            ndimf  = obj.dim.ndimField;
            nnode  = obj.dim.nnode;
            Ae     = aElem;
            A = sparse(ndofs,ndofs);
            for i = 1:nnode*ndimf
                dofsI = dofConnec(:,i);
                for j = 1:nnode*ndimf
                    dofsJ = dofConnec(:,j);
                    a = squeeze(Ae(i,j,:));
                    Aadd = obj.computeAaddBySparse(a, dofsI, dofsJ);
                    Aadd = obj.computeAaddByAccumarray(a, dofsI, dofsJ);
                    A = A + Aadd;
                end
            end
        end

        function A = assembleMatrixViaIndices3(obj, aElem)
            disp('(--calculem DOF connec')
            tic
            dofConnec = obj.computeDofConnectivity()';
            toc
            disp('--)')
%             dofConnec = obj.computeDofConnectivity()';
            ndimf  = obj.dim.ndimField;
            nnode  = obj.dim.nnode;
            Ae     = aElem;
%             A = sparse(ndofs,ndofs);
            res = [];
            for i = 1:nnode*ndimf
                dofsI = dofConnec(:,i);
                for j = 1:nnode*ndimf
                    dofsJ = dofConnec(:,j);
                    a = squeeze(Ae(i,j,:));
                    res = [res; dofsI, dofsJ, a];
                end
            end
            A = sparse(res(:,1), res(:,2), res(:,3));
        end

        function A = assembleMatrixViaIndices(obj, aElem)
            dofConnec = obj.computeDofConnectivity()';
            ndimf  = obj.dim.ndimField;
            nnode  = obj.dim.nnode;
            Ae     = aElem;
            ndofEl = nnode*ndimf;
            ndof = obj.dim.ndof;
            nelem  = size(dofConnec, 1); % obj.dim.nelem does NOT work
%             A = sparse(ndofs,ndofs);
            res = zeros(ndofEl^2 * nelem, 3);
            strt = 1;
            fnsh = nelem;
            for i = 1:ndofEl
                dofsI = dofConnec(:,i);
                for j = 1:ndofEl
                    dofsJ = dofConnec(:,j);
                    a = squeeze(Ae(i,j,:));
                    matRes = [dofsI, dofsJ, a];
                    res(strt:fnsh,:) = matRes;
                    strt = strt + nelem;
                    fnsh = fnsh + nelem;
                end
            end
            A = sparse(res(:,1), res(:,2), res(:,3), ndof, ndof);
        end
        
        function dofConnec = computeDofConnectivity(obj)
            connec = obj.globalConnec;
            ndimf  = obj.dim.ndimField;
            nnode  = obj.dim.nnode;
            dofsElem  = zeros(nnode*ndimf,size(connec,1));
            for inode = 1:nnode
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
            ndimf = obj.dim.ndimField;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

        function Cadd = computeAaddBySparse(obj,a, dofsI, dofsJ)
           d = obj.dim;
           ndofs = d.ndof;
           Cadd = sparse(dofsI,dofsJ,a,ndofs,ndofs);
        end

        function Cadd = computeAaddByAccumarray(obj,a, dofsI, dofsJ)
           d = obj.dim;
           ndofs = d.ndof;
           index = [dofsI, dofsJ];
           Cadd = accumarray(index,a,[ndofs ndofs],[],[],true);
        end

        %% Vector assembly

        function V = assembleVector(obj, F, dofsInElem)
            dofsInElem = obj.computeDofConnectivity();
            ndofPerElem = obj.dim.ndofPerElement;
            ndof        = obj.dim.ndof;
            V = zeros(ndof,1);
            for iDof = 1:ndofPerElem
                dofs = dofsInElem(iDof,:);
                c = F(iDof,:);
                Fadd = obj.computeAddVectorBySparse(dofs, c);
                % Fadd = obj.computeAddVectorByAccumarray(dofs, c);
                V = V + Fadd;
            end
        end

        function Vadd = computeAddVectorBySparse(obj,dofs, c)
           ndof = obj.dim.ndof;
           Vadd = sparse(dofs,1,c',ndof,1);
        end

        function Vadd = computeAddVectorByAccumarray(obj,dofs,c)
           ndof = obj.dim.ndof;
           Vadd = accumarray(dofs',c',[ndof 1]);
        end
        
        %% B matrix assembly
        function Bt = assembleBMatrix(obj, Bfull)
            d = obj.dim;
            ntot  = d.nt;
            ndofGlob = d.ndof;
            Bt = sparse(ntot,ndofGlob);
            for idof = 1:d.ndofPerElement
                dofs  = obj.computeGlobalDofs(idof);
                Bidof = Bfull(:,idof);
                Bdof = obj.computeBdofBySparse(Bidof,dofs);
                %                Bdof = obj.computeBdofByAccumarray(Bidof,dofs);
                Bt = Bt + Bdof;
            end
        end

        function dofs = computeGlobalDofs(obj,idof)
            d = obj.dim;
            ngaus = d.ngaus;
            nstre = d.nstre;
            gDofs = obj.transformLocal2Global(idof);
            dofs = repmat(gDofs',ngaus*nstre,1);
            dofs = dofs(:);
        end

        function gDofs = transformLocal2Global(obj,iDof)
            d     = obj.dim;
            ndimf = d.ndimField;
            nnode = d.nnode;
            nodes        = obj.globalConnec;
            nodesInElem  = reshape(repmat(1:nnode,ndimf,1),1,[]);
            dofs         = repmat(1:ndimf,1,nnode);
            inode        = nodesInElem(iDof);
            iunkn        = dofs(iDof);
            nodeI        = nodes(:,inode);
            gDofs   = ndimf*(nodeI-1) + iunkn;
        end

        function Bdof = computeBdofBySparse(obj,Bidof,dofs)
            d = obj.dim;
            ntot  = d.nt;
            ndofGlob = d.ndof;
            Bdof = sparse(1:ntot,dofs,Bidof,ntot,ndofGlob);
        end

        function Bdof = computeBdofByAccumarray(obj,Bidof,dofs)
            d = obj.dim;
            ntot  = d.nt;
            ndof = d.ndof;
            posI = 1:ntot;
            index = [posI', dofs];
            %            Bdof = accumarray(dofs,Bidof,[ntot 1]);
            Bdof = accumarray(index,Bidof,[ntot ndof],[],[],true);
        end

        %% C matrix assembly

       function CmatTot = assembleCMatrix(obj, Cmat, dvol)
           nstre = obj.dim.nstre;
           ngaus = obj.dim.ngaus;
           ntot  = obj.dim.nt;
           CmatTot = sparse(ntot,ntot);
           for istre = 1:nstre
               for jstre = 1:nstre
                   for igaus = 1:ngaus
                       posI = (istre)+(nstre)*(igaus-1) : ngaus*nstre : ntot;
                       posJ = (jstre)+(nstre)*(igaus-1) : ngaus*nstre : ntot;
                       Ci = Cmat(istre,jstre,:,igaus);
                       Ct = squeeze(Ci).*dvol(:,igaus);
                       Cadd = obj.computeCaddBySparse(Ct, posI, posJ);
                       CmatTot = CmatTot + Cadd;
                   end
               end
           end
       end

       function Cadd = computeCaddBySparse(obj,Ct, posI, posJ)
           d = obj.dim;
           ntot  = d.nt;
           Cadd = sparse(posI,posJ,Ct,ntot,ntot);
       end

       function Cadd = computeCaddByAccumarray(obj,Ct, posI, posJ)
           d = obj.dim;
           ntot  = d.nt;
           index = [posI', posJ'];
           Cadd = accumarray(index,Ct,[ntot ntot],[],[],true);
      end


    end

end