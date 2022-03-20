classdef Assembler < handle

    properties (Access = private)
        dim
        globalConnec
        dofsInElem
    end

    methods (Access = public)
        
        function obj = Assembler(cParams)
            obj.init(cParams);
        end

        function A = assemble(obj, aElem)
            A = obj.assembleMatrix(aElem);
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
            obj.dofsInElem   = cParams.dofsInElem;
        end

        function A = assembleMatrix(obj, aElem)
            dofConnec = obj.dofsInElem';
            ndofs  = obj.dim.ndof;
            ndimf  = obj.dim.ndimField;
            nnode  = obj.dim.nnode;
            Ae     = aElem;
            A = sparse(ndofs,ndofs);
            assemblenum1 = size(aElem,1);
            assemblenum2 = size(aElem,2);
            for i = 1:nnode*ndimf
                dofsI = dofConnec(:,i);
                for j = 1:nnode*ndimf
                    dofsJ = dofConnec(:,j);
                    a = squeeze(Ae(i,j,:));
                    Aadd = obj.computeAaddBySparse(a, dofsI, dofsJ);
%                     Aadd = obj.computeAaddByAccumarray(a, dofsI, dofsJ);
                    A = A + Aadd;
                end
            end
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