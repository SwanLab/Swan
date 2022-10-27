classdef BMatrixComputer < handle

    properties (Access = private)
        dim
        nvoigt
        geometry
    end

    methods (Access = public)

        function obj = BMatrixComputer(cParams)
            obj.init(cParams);
        end

        function B = compute(obj,igaus)
            ndimf = obj.dim.ndimf;
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
        end

        function B = computeBin2D(obj,igaus)
            d = obj.dim;
            dNdx = obj.geometry.dNdx;
            nstre = obj.nvoigt;
            ndimf = d.ndimf;
            nnode = size(dNdx,2);
            ndofE = nnode*ndimf;
            nelem = size(dNdx,3);
            B = zeros(nstre,ndofE,nelem);
            for inode = 1:nnode
                j = ndimf*(inode-1)+1;
                B(1,j,:)   = dNdx(1,inode,:,igaus);
                B(2,j+1,:) = dNdx(2,inode,:,igaus);
                B(3,j,:)   = dNdx(2,inode,:,igaus);
                B(3,j+1,:) = dNdx(1,inode,:,igaus);
            end
        end

        function B = computeBin3D(obj,igaus)
            d    = obj.dim;
            dNdx = obj.geometry.dNdx;
            nNode = size(dNdx,2);
            nElem = size(dNdx,3);
            B = zeros(obj.nvoigt,nNode,nElem);
            for inode = 1:nNode
                j = d.ndimf*(inode-1)+1;
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
            dNdx  = obj.geometry.dNdx(:,:,:,igaus);
            d     = obj.dim;
            nNode = size(dNdx,2);
            nElem = size(obj.geometry.dNdx,3);
            nDofs = d.ndimf*nNode;
            B = zeros(obj.nvoigt,nDofs,nElem);
            for inode = 1:nNode
                j = d.ndimf*(inode-1) + 1;
                B(1,j,:) = dNdx(1,inode,:);
                B(2,j,:) = dNdx(2,inode,:);
            end
        end

    end
end

