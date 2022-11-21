classdef LHSintegrator_StiffnessBeam < LHSintegrator

    properties (Access = private)
        coord
        connec
        material
    end

    methods (Access = public)

        function obj = LHSintegrator_StiffnessBeam(cParams)
            obj.coord    = cParams.coord;
            obj.connec   = cParams.connec;
            obj.material = cParams.material;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrixField(lhs);
        end

    end


    methods (Access = private)

        function lhs = computeElementalLHS(obj)
            le = obj.calculateBarLength();
            [A, Iz] = obj.getBarSectionData();
            Ke = obj.computeLocalK(le, A);
            Re = obj.computeRotationMatrix(le);
            ReT   = permute(Re,[2 1 3]);
            ReTKe  = pagemtimes(ReT,Ke);
            lhs    = pagemtimes(ReTKe, Re);
        end

        function LHS = assembleMatrixField(obj, Ae)
            nNods = size(obj.coord, 1);
            nNodE = size(obj.connec,2);
            nDim  = 3;
            nDofs = nNods*nDim;
            dofConn = obj.computeDofConnectivity();
            Kg = zeros(nNods*nDim, nNods*nDim);
            for i = 1:nNodE*nDim
                for j = 1:nNodE*nDim
                    a = squeeze(Ae(i,j,:));
                    Kg = Kg + sparse(dofConn(i,:),dofConn(j,:),a,nDofs,nDofs);
                end
            end
            LHS = Kg;
        end

        function le = calculateBarLength(obj)
            [dX, dY, dZ] = obj.calculateBarDeltas();
            le = sqrt(dX.^2 + dY.^2 + dZ.^2);
        end

        function [dX, dY, dZ] = calculateBarDeltas(obj)
            nBars = size(obj.connec,1);
            nod1 = obj.connec(:,1);
            nod2 = obj.connec(:,2);
            x1 = obj.coord(nod1, 1);
            x2 = obj.coord(nod2, 1);
            y1 = obj.coord(nod1, 2);
            y2 = obj.coord(nod2, 2);
            z1 = obj.coord(nod1, 2);
            z2 = obj.coord(nod2, 2);
            dX = reshape(x2 - x1, [1,1,nBars]);
            dY = reshape(y2 - y1, [1,1,nBars]);
            dZ = reshape(z2 - z1, [1,1,nBars]);
        end

        function [A, Iz] = getBarSectionData(obj)
            nBars = size(obj.connec,1);
            A  = ones(nBars, 1);
            Iz = ones(nBars, 1);
            A  = reshape(A, [1,1, nBars]);
            Iz = reshape(Iz, [1,1, nBars]);
        end

        function ke = computeLocalK(obj, le, A)
            E = obj.material.E;
            nBars = size(obj.connec,1);
            mat = ones(2,2, nBars);
            mat(1,2, :) = -1;
            mat(2,1, :) = -1;
            ke = E*A./le.* mat;
        end

        function re = computeRotationMatrix(obj, le)
            nBars = size(obj.connec,1);
            [dX, dY, dZ] = obj.calculateBarDeltas();
            coef = 1./le;
            re = zeros(2,6, nBars);
            re(1,1,:) = coef.*dX;
            re(1,2,:) = coef.*dY;
            re(1,3,:) = coef.*dZ;
            re(2,4,:) = coef.*dX;
            re(2,5,:) = coef.*dY;
            re(2,6,:) = coef.*dZ;
        end

        function dofConnec = computeDofConnectivity(obj)
            conne = obj.connec;
            nDimf = size(obj.coord,2);
            nNodeEl = size(conne, 2);
            nDofsEl = nNodeEl * nDimf;
            dofsElem  = zeros(nDofsEl,size(conne,1));
            for inode = 1:nNodeEl
                for iunkn = 1:nDimf
                    idofElem   = obj.nod2dof(inode,iunkn);
                    globalNode = conne(:,inode);
                    idofGlobal = obj.nod2dof(globalNode,iunkn);
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            dofConnec = dofsElem';
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = size(obj.coord,2);
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

    end

end