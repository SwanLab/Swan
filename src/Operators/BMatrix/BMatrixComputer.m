classdef BMatrixComputer < handle

    properties (Access = private)
        fun
        dNdx
        nVoigt
    end

    methods (Access = public)

        function obj = BMatrixComputer(cParams)
            obj.init(cParams);
        end

        function B = compute(obj,igaus)
            ndimf = obj.fun.ndimf;
            switch ndimf
                case 1
                    obj.nVoigt = 2;
                    B = obj.computeBin1D(igaus);
                case 2
                    obj.nVoigt = 3;
                    B = obj.computeBin2D(igaus);
                case 3
                    obj.nVoigt = 6;
                    B = obj.computeBin3D(igaus);
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.fun  = cParams.fun;
            obj.dNdx = cParams.dNdx;
        end

        function B = computeBin2D(obj,iGaus)
            deriv = obj.dNdx;
            nStre = obj.nVoigt;
            nDimf = obj.fun.ndimf;
            nNodE = size(deriv,2);
            nDofE = nNodE*nDimf;
            nElem = size(deriv,4);
            B = zeros(nStre,nDofE,nElem);
            for iNode = 1:nNodE
                j = nDimf*(iNode-1)+1;
                B(1,j,:)   = deriv(1,iNode,iGaus,:);
                B(2,j+1,:) = deriv(2,iNode,iGaus,:);
                B(3,j,:)   = deriv(2,iNode,iGaus,:);
                B(3,j+1,:) = deriv(1,iNode,iGaus,:);
            end
        end

        function B = computeBin3D(obj,iGaus)
            deriv = obj.dNdx;
            nNode = size(deriv, 2);
            nElem = size(deriv, 4);
            j = reshape((0:nNode-1) * obj.fun.ndimf, 1, nNode, 1);
            d1 = permute(deriv(1, :, iGaus, :), [1, 2, 4, 3]);
            d2 = permute(deriv(2, :, iGaus, :), [1, 2, 4, 3]);
            d3 = permute(deriv(3, :, iGaus, :), [1, 2, 4, 3]);
            B = zeros(obj.nVoigt, obj.fun.ndimf * nNode, nElem);
            B(1, j + 1, :) = d1;
            B(2, j + 2, :) = d2;
            B(3, j + 3, :) = d3;
            B(4, j + 1, :) = d2;
            B(4, j + 2, :) = d1;
            B(5, j + 1, :) = d3;
            B(5, j + 3, :) = d1;
            B(6, j + 2, :) = d3;
            B(6, j + 3, :) = d2;
        end

        function [B] = computeBin1D(obj, iGaus)
            deriv  = obj.dNdx;
            nDimf = obj.fun.ndimf;
            nNode = size(deriv,2);
            nElem = size(obj.dNdx,4);
            nDofs = nDimf*nNode;
            B = zeros(obj.nVoigt,nDofs,nElem);
            for inode = 1:nNode
                j = nDimf*(inode-1) + 1;
                B(1,j,:) = deriv(1,inode,iGaus,:);
                B(2,j,:) = deriv(2,inode,iGaus,:);
            end
        end

    end
end

