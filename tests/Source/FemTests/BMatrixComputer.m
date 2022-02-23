classdef BMatrixComputer < handle

    properties (Access = private)
        dim
        geometry
        globalConnec
    end

    methods (Access = public)

        function obj = BMatrixComputer(cParams)
            obj.init(cParams);
        end

        function Btot = compute(obj)
            Bmatrix  = obj.computeBinMatrixForm();
            Btot     = obj.assembleMatrix(Bmatrix);
        end

        function B = computeBmat(obj,igaus)
            ndim = obj.dim.ndim;
            switch ndim
                case 2
                    B = obj.computeBin2D(igaus);
                case 3
                    B = obj.computeBin3D(igaus);
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.dim          = cParams.dim;
            obj.geometry     = cParams.geometry;
            obj.globalConnec = cParams.globalConnec;
        end

        function B = computeBin2D(obj,igaus)
            d = obj.dim;
            nstre          = d.nstre;
            nnode          = d.nnode;
            nelem          = d.nelem;
            nunkn          = d.nunkn;
            ndofPerElement = d.ndofPerElement;
            dNdx = obj.geometry.dNdx;
            B = zeros(nstre,ndofPerElement,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)   = dNdx(1,i,:,igaus);
                B(2,j+1,:) = dNdx(2,i,:,igaus);
                B(3,j,:)   = dNdx(2,i,:,igaus);
                B(3,j+1,:) = dNdx(1,i,:,igaus);
            end
        end

        function B = computeBin3D(obj,igaus)
            d    = obj.dim;
            dNdx = obj.geometry.dNdx;
            B = zeros(d.nstre,d.ndofPerElement,d.nelem);
            for inode=1:d.nnode
                j = d.nunkn*(inode-1)+1;
                % associated to normal strains
                B(1,j,:)   = dNdx(1,inode,:,igaus);
                B(2,j+1,:) = dNdx(2,inode,:,igaus);
                B(3,j+2,:) = dNdx(3,inode,:,igaus);
                % associated to shear strain, gamma12
                B(4,j,:)   = dNdx(2,inode,:,igaus);
                B(4,j+1,:) = dNdx(1,inode,:,igaus);
                % associated to shear strain, gamma13
                B(5,j,:)   = dNdx(3,inode,:,igaus);
                B(5,j+2,:) = dNdx(1,inode,:,igaus);
                % associated to shear strain, gamma23
                B(6,j+1,:) = dNdx(3,inode,:,igaus);
                B(6,j+2,:) = dNdx(2,inode,:,igaus);
            end
        end

        function Bmatrix = computeBinMatrixForm(obj)
            d  = obj.dim;
            nB = d.nstre*d.ngaus*d.nelem;
            Bmatrix = zeros(nB,d.ndofPerElement);
            for igaus = 1:d.ngaus
                Bgaus = obj.computeBmatrix(igaus);
                index = obj.computeGlobalIndex(igaus);
                Bmatrix(index,:) = Bgaus;
            end
        end

        function B = computeBmatrix(obj,igaus)
            d  = obj.dim;
            Bmat = obj.computeBmat(igaus);
            Bper = permute(Bmat,[1 3 2]);
            B    = reshape(Bper,d.nelem*d.nstre,d.ndofPerElement);
        end

        function index = computeGlobalIndex(obj,igaus)
            d = obj.dim;
            uIndex = obj.computeUnitaryIndex(igaus);
            index = repmat(uIndex,d.nelem,1);
        end

        function index = computeUnitaryIndex(obj,igaus)
            d = obj.dim;
            nGaus = d.ngaus;
            nStre = d.nstre;
            index = false(nGaus*nStre,1);
            pos =  nStre*(igaus-1) + (1:nStre);
            index(pos) = true;
        end

        function Bt = assembleMatrix(obj, Bfull)
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
            nunkn = d.nunkn;
            nodes        = obj.globalConnec;
            nodesInElem  = reshape(repmat(1:d.nnode,d.nunkn,1),1,[]);
            dofs         = repmat(1:d.nunkn,1,d.nnode);
            inode        = nodesInElem(iDof);
            iunkn        = dofs(iDof);
            nodeI        = nodes(:,inode);
            gDofs   = nunkn*(nodeI-1) + iunkn;
        end

    end
end

