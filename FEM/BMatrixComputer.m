classdef BMatrixComputer < handle

    properties (Access = private)
        dim
        nvoigt
        geometry
        globalConnec
    end

    methods (Access = public)

        function obj = BMatrixComputer(cParams)
            obj.init(cParams);
        end

        function Btot = compute(obj)
            Bmatrix = obj.computeBinMatrixForm();
            Btot    = obj.assembleMatrix(Bmatrix);
        end

        function B = computeBmat(obj,igaus)
            ndimf = obj.dim.ndimField;
            switch ndimf
                case 1
                    obj.nvoigt = 2;
                    B = obj.computeBin1D(igaus);
                case 2
                    obj.nvoigt = 3;
                    B = obj.computeBin2D(igaus);
                case 3
                    obj.nvoigt = 6;
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
            dNdx = obj.geometry.dNdx;
            nstre = obj.nvoigt;
            nnode = d.nnode;
            ndimf = d.ndimField;
            ndofE = d.ndofPerElement;
            nelem = size(dNdx,3);
            B = zeros(nstre,ndofE,nelem);
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
            nelem = size(dNdx,3);
            B = zeros(obj.nvoigt,d.ndofPerElement,nelem);
            for inode = 1:d.nnode
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
            d     = obj.dim;
            nelem = size(obj.geometry.dNdx,3);
            dNdx  = obj.geometry.dNdx(:,:,:,igaus);
            B = zeros(obj.nvoigt,d.nnode*d.ndimField,nelem);
            for inode = 1:d.nnode
                j = d.ndimField*(inode-1) + 1;
                B(1,j,:) = dNdx(1,inode,:);
                B(2,j,:) = dNdx(2,inode,:);
            end
        end

        function Bmatrix = computeBinMatrixForm(obj)
            nelem = size(obj.geometry.dNdx,3);
            ngaus = size(obj.geometry.dNdx,4);
            ndofE = obj.dim.ndofPerElement;
            nB = obj.nvoigt*ngaus*nelem;
            Bmatrix = zeros(nB,ndofE);
            for igaus = 1:ngaus
                Bgaus = obj.computeBmatrix(igaus);
                index = obj.computeGlobalIndex(igaus);
                Bmatrix(index,:) = Bgaus;
            end
        end

        function B = computeBmatrix(obj,igaus)
            nelem = size(obj.geometry.dNdx,3);
            ndofE = obj.dim.ndofPerElement;
            Bmat = obj.computeBmat(igaus);
            Bper = permute(Bmat,[1 3 2]);
            B    = reshape(Bper,nelem*obj.nvoigt,ndofE);
        end

        function index = computeGlobalIndex(obj,igaus)
            nelem = size(obj.geometry.dNdx,3);
            uIndex = obj.computeUnitaryIndex(igaus);
            index = repmat(uIndex,nelem,1);
        end

        function index = computeUnitaryIndex(obj,igaus)
            nGaus = size(obj.geometry.dNdx,4);
            nstre = obj.nvoigt;
            index = false(nGaus*nstre,1);
            pos =  nstre*(igaus-1) + (1:nstre);
            index(pos) = true;
        end

        function Bt = assembleMatrix(obj, Bfull)
            s.dim = obj.dim;
            s.globalConnec = obj.globalConnec;
            dNdx = obj.geometry.dNdx;
            d.nelem = size(dNdx,3);
            d.ngaus = size(dNdx,4);
            assembler = Assembler(s);
            Bt = assembler.assembleB(Bfull, d);
        end

    end
end

