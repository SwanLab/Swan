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

        function [Btot, Balt] = compute(obj)
            Bmatrix = obj.computeBinMatrixForm();
            Btot    = obj.assembleMatrix(Bmatrix);
            Balt = 0;
        end

        function B = computeBmat(obj,igaus)
            ndim = obj.dim.ndim;
            switch ndim
                case 1
                    B = obj.computeBin1D(igaus);
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
            ndimf          = d.ndimField;
            ndofPerElement = d.ndofPerElement;
            dNdx = obj.geometry.dNdx;
            B = zeros(nstre,ndofPerElement,nelem);
            for i = 1:nnode
                j = ndimf*(i-1)+1;
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
                j = d.ndimField*(inode-1)+1;
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

        function [B] = computeBin1D(obj, igaus)
            d    = obj.dim;
            dNdx = obj.geometry.dNdx(:,:,:,igaus);
            B = zeros(2,d.nnode*d.ndimField,d.nelem);
            for inode = 1:d.nnode
                j = d.ndimField*(inode-1) + 1;
                B(1,j,:) = dNdx(1,inode,:);
                B(2,j,:) = dNdx(2,inode,:);
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
            s.dim = obj.dim;
            s.globalConnec = obj.globalConnec;
            assembler = Assembler(s);
            Bt = assembler.assembleB(Bfull);
        end


        function Bmatrix = computeAlternateB(obj)
            d  = obj.dim;
            nB = d.nstre*d.ngaus*d.nelem;
            Bmatrix = zeros(nB,d.ndofPerElement);
            for igaus = 1:d.ngaus
                Bmat = obj.computeBmat(igaus);
            end
        end

    end
end

