classdef LHSintegrator_Stiffness3DBeam < LHSintegrator

    properties (Access = private)
        coord
        connec
        material
    end

    methods (Access = public)

        function obj = LHSintegrator_Stiffness3DBeam(cParams)
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
            [A, Iz] = obj.interp.computeSectionAreaAndInertia();
            Ke = obj.computeLocalK(le, A, Iz);
            Re = obj.computeRotationMatrix(le);
            ReT   = permute(Re,[2 1 3]);
            ReTKe  = pagemtimes(ReT,Ke);
            lhs    = pagemtimes(ReTKe, Re);
        end

        function lhs = computeElementalLHSDerivative(obj)
            le = obj.calculateBarLength();
            [dA, dIz] = obj.interp.computeSectionAreaAndInertiaDerivative();
            Ke = obj.computeLocalKDerivative(le, dA, dIz);
            Re = obj.computeRotationMatrix(le);
            ReT   = permute(Re,[2 1 3]);
            ReTKe  = pagemtimes(ReT,Ke);
            lhs    = pagemtimes(ReTKe, Re);
        end 

        function LHS = assembleMatrixField(obj, Ae)
            nNods = size(obj.coord, 1);
            nNodE = size(obj.connec,2);
            nDim  = 6; %3
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

        function Ke = computeLocalK(obj, L, A, Iz)
            E = obj.material.E;
            G = obj.material.G;
            % J = St Venant torsional stiffness
            Iy = Iz;
            J = Iz/2; % For ISCSO sections
            nBars = size(obj.connec,1);
            Ke = zeros(12,12, nBars);
            k1 = E.*A./L;
            k2 = G.*J./L;

            Ke(1,1,:) = +k1;
            Ke(1,7,:) = -k1;
            Ke(7,1,:) = -k1;
            Ke(7,7,:) = +k1;

            Ke(2,2,:) = +12.*E.*Iz./L.^3;
            Ke(2,8,:) = -12.*E.*Iz./L.^3;
            Ke(8,2,:) = -12.*E.*Iz./L.^3;
            Ke(8,8,:) = +12.*E.*Iz./L.^3;

            Ke(2,6,:)  = +6.*E.*Iz./L.^2;
            Ke(2,12,:) = +6.*E.*Iz./L.^2;
            Ke(8,6,:)  = -6.*E.*Iz./L.^2;
            Ke(8,12,:) = -6.*E.*Iz./L.^2;

            Ke(3,3,:) = +12.*E.*Iy./L.^3;
            Ke(3,9,:) = -12.*E.*Iy./L.^3;
            Ke(9,3,:) = -12.*E.*Iy./L.^3;
            Ke(9,9,:) = +12.*E.*Iy./L.^3;

            Ke(3,5,:)  = -6.*E.*Iy./L.^2;
            Ke(3,11,:) = -6.*E.*Iy./L.^2;
            Ke(9,5,:)  = +6.*E.*Iy./L.^2;
            Ke(9,11,:) = +6.*E.*Iy./L.^2;

            Ke(4,4,:)   = +k2;
            Ke(4,10,:)  = -k2;
            Ke(10,4,:)  = -k2;
            Ke(10,10,:) = +k2;

            Ke(5,3,:)  = -6.*E.*Iy./L.^2;
            Ke(5,9,:)  = +6.*E.*Iy./L.^2;
            Ke(11,3,:) = -6.*E.*Iy./L.^2;
            Ke(11,9,:) = +6.*E.*Iy./L.^2;

            Ke(5,5,:)   = +4.*E.*Iy./L;
            Ke(5,11,:)  = +2.*E.*Iy./L;
            Ke(11,5,:)  = +2.*E.*Iy./L;
            Ke(11,11,:) = +4.*E.*Iy./L;

            Ke(6,2,:)  = +6.*E.*Iz./L.^2;
            Ke(6,8,:)  = -6.*E.*Iz./L.^2;
            Ke(12,2,:) = +6.*E.*Iz./L.^2;
            Ke(12,8,:) = -6.*E.*Iz./L.^2;

            Ke(6,6,:)   = +4.*E.*Iz./L;
            Ke(6,12,:)  = +2.*E.*Iz./L;
            Ke(12,6,:)  = +2.*E.*Iz./L;
            Ke(12,12,:) = +4.*E.*Iz./L;
        end

        function Re = computeRotationMatrix(obj, L)
            nBars = size(obj.connec,1);
            [dX, dY, dZ] = obj.calculateBarDeltas();
            ReT = zeros(3,3, nBars);
            Re0 = zeros(3,3, nBars);

            % https://www.solid.lth.se/fileadmin/hallfasthetslara/utbildning/kurser/FHL064_FEM/calfem34.pdf
            % Page 26
            % nxX specifies the cosine of the angle between the x axis and
            % X axis, and so on.

            % https://www.brown.edu/Departments/Engineering/Courses/En2340/Projects/Projects_2015/Wenqiang_Fan.pdf
            % Page 5

            T11 = dX./L; % these are the easy ones
            T12 = dY./L; % coz the cos are with X
            T13 = dZ./L; % and we already know X (beam direction)

            % we need to define a vector k that is perpendicular to X

            xVec   = [dX, dY, dZ];
            jkVecs = null(xVec);
            kVec   = jkVecs(:,2); 
            k1 = kVec(1);
            k2 = kVec(2);
            k3 = kVec(3);

            A = sqrt( (T12*k3 - T13*k2).^2 + (T13*k1 - T11*k3).^2 + ...
                (T11*k2 - T12*k1).^2);

            T21 = -(T12*k3 - T13*k2)./A;
            T22 = -(T13*k1 - T11*k3)./A;
            T23 = -(T11*k2 - T12*k1)./A;

            T31 = +(T12*T23 - T13*T22)./B;
            T32 = -(T13*T21 - T11*T23)./B;
            T33 = -(T11*T22 - T12*T21)./B;

            ReT(1,1,:) = T11; % nxX;
            ReT(1,2,:) = T12; % nyX;
            ReT(1,3,:) = T13; % nzX;

            ReT(2,1,:) = T21; % nxY;
            ReT(2,2,:) = T22; % nyY;
            ReT(2,3,:) = T23; % nzY;

            ReT(3,1,:) = T31; % nxZ;
            ReT(3,2,:) = T32; % nyZ;
            ReT(3,3,:) = T33; % nzZ;
            Re = [ReT, Re0, Re0, Re0;
                  Re0, ReT, Re0, Re0;
                  Re0, Re0, ReT, Re0;
                  Re0, Re0, Re0, ReT;
                 ];
        end

        function dofConnec = computeDofConnectivity(obj)
            conne = obj.connec;
            nDimf = 6;
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
            ndimf = 6;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
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

    end

end