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

        function LHS = computeDerivative(obj)
            lhs = obj.computeElementalLHSDerivative();
            LHS = obj.assembleMatrixField(lhs);
        end

    end


    methods (Access = private)

        function lhs = computeElementalLHS(obj)
            le = obj.calculateBarLength();
            [A, Iz] = obj.getBarSectionData();
            Ke = obj.computeLocalK(le, A, Iz);
            Re = obj.computeRotationMatrix(le);
            ReT   = permute(Re,[2 1 3]);
            ReTKe  = pagemtimes(ReT,Ke);
            lhs    = pagemtimes(ReTKe, Re);
        end

        function lhs = computeElementalLHSDerivative(obj)
            le = obj.calculateBarLength();
            [dA, dIz] = obj.getBarDerivativeSectionData();
            Ke = obj.computeLocalKDerivative(le, A, Iz);
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
            for i = 1:size(Ae,1)
                for j = 1:size(Ae,1)
                    a = squeeze(Ae(i,j,:));
                    Kg = Kg + sparse(dofConn(i,:),dofConn(j,:),a,nDofs,nDofs);
                end
            end
            LHS = Kg;
        end

        function le = calculateBarLength(obj)
            [dX, dY] = obj.calculateBarDeltas();
%             le = sqrt(dX.^2 + dY.^2 + dZ.^2);
            le = sqrt(dX.^2 + dY.^2);
        end

        function [dX, dY] = calculateBarDeltas(obj)
            nBars = size(obj.connec,1);
            nod1 = obj.connec(:,1);
            nod2 = obj.connec(:,2);
            x1 = obj.coord(nod1, 1);
            x2 = obj.coord(nod2, 1);
            y1 = obj.coord(nod1, 2);
            y2 = obj.coord(nod2, 2);
%             z1 = obj.coord(nod1, 2);
%             z2 = obj.coord(nod2, 2);
            dX = reshape(x2 - x1, [1,1,nBars]);
            dY = reshape(y2 - y1, [1,1,nBars]);
%             dZ = reshape(z2 - z1, [1,1,nBars]);
        end

        function [A, Iz] = getBarSectionData(obj)
            nBars = size(obj.connec,1);
            A  = ones(nBars, 1);
            Iz = ones(nBars, 1);
            A  = reshape(A, [1,1, nBars]);
            Iz = reshape(Iz, [1,1, nBars]);
        end

        function Ke = computeLocalK(obj, le, A, Iz)
            E = obj.material.E;
            nBars = size(obj.connec,1);
            Ke = zeros(6,6, nBars);
            c1 = Iz.*E./le.^3;
            c2 = A.*E./le;
            Ke(1,1,:) = c2;
            Ke(1,4,:) = -c2;
            Ke(2,2,:) = c1*12;
            Ke(2,3,:) = c1.*6.*le;
            Ke(2,5,:) = c1*(-12);
            Ke(2,6,:) = c1.*6.*le;
            Ke(3,2,:) = c1.*6.*le;
            Ke(3,3,:) = c1.*4.*le.^2;
            Ke(3,5,:) = c1.*(-6*le);
            Ke(3,6,:) = c1.* 2.*le.^2;
            Ke(4,1,:) = -c2;
            Ke(4,4,:) = c2;
            Ke(5,2,:) = -12*c1;
            Ke(5,3,:) = -6*le.*c1;
            Ke(5,5,:) = 12*c1;
            Ke(5,6,:) = -6*le.*c1;
            Ke(6,2,:) = 6*le.*c1;
            Ke(6,3,:) = 2*le.^2.*c1;
            Ke(6,5,:) = -6.*le.*c1;
            Ke(6,6,:) = 4.*le.^2.*c1;

        end

        function dKe = computeLocalKDerivative(obj, le, dA, dIz)
            dAdr = dA(1);
            dAdt = dA(2);
            dIdr = dIz(1);
            dIdt = dIz(2);
            dKedr = computeLocalK(obj, le, dAdr, dIdr);
            dKedt = computeLocalK(obj, le, dAdt, dIdt);
            dKe(:,:,:,1) = dKedr;
            dKe(:,:,:,2) = dKedt;
            % suposo que els termes creuats no calen?
        end

        function Re = computeRotationMatrix(obj, le)
            nBars = size(obj.connec,1);
            [dX, dY] = obj.calculateBarDeltas();
            Re = zeros(6,6, nBars);
            coef = 1./le;
            Re(1,1,:) =   coef.*dX;
            Re(1,2,:) =   coef.*dY;
            Re(2,1,:) = - coef.*dY;
            Re(2,2,:) =   coef.*dX;
            Re(3,3,:) =   coef.*le;
            Re(4,4,:) =   coef.*dX;
            Re(4,5,:) =   coef.*dY;
            Re(5,4,:) = - coef.*dY;
            Re(5,5,:) =   coef.*dX;
            Re(6,6,:) =   coef.*le;
        end

        function dofConnec = computeDofConnectivity(obj)
            conne = obj.connec;
            nDimf = 3;
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
            ndimf = 3;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

    end

end