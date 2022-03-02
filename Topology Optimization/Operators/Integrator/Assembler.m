classdef Assembler < handle

    properties (Access = private)
        dim
        globalConnec
    end

    methods (Access = public)
        
        function obj = Assembler(cParams)
            obj.init(cParams);
        end

        function A = assemble(obj, aElem)
            A = obj.assembleMatrix(aElem);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.dim          = cParams.dim;
            obj.globalConnec = cParams.globalConnec;
        end

        function A = assembleMatrix(obj, aElem)
            connec    = obj.globalConnec;
            dofConnec = obj.computeDOFconnec();
            ndofs  = obj.dim.ndof;
            ndimf  = obj.dim.ndimField;
%             nnode  = obj.dim.nnode;
            nnode  = size(connec,2); % pending review why TopOptTests take incorrect nnode
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

        function dof_elem = computeDOFconnec(obj)
            connec = obj.globalConnec;
            ndimf  = obj.dim.ndimField;
            nnode  = size(connec,2);
            dof_elem  = zeros(nnode*ndimf,size(connec,1));
            for inode = 1:nnode
                for iunkn = 1:ndimf
                    idof_elem = obj.nod2dof(inode,iunkn);
                    global_node = connec(:,inode);
                    idof_global = obj.nod2dof(global_node,iunkn);
                    dof_elem(idof_elem,:) = idof_global;
                end
            end
            dof_elem = dof_elem';
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimField;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

    end

end