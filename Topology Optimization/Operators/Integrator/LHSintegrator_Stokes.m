classdef LHSintegrator_Stokes < handle %LHSintegrator

    properties (GetAccess = public, SetAccess = private)
        Melem
        Kelem
        Delem, D
    end

    properties (Access = private)
        dt
        mesh
        velocityField
        pressureField
    end

    methods (Access = public)

        function obj = LHSintegrator_Stokes(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj)
            velLHS = obj.computeVelocityLHS();
            D      = obj.computeDmatrix();
            prsLHS = obj.computePressureLHS();
            LHS = [velLHS, D; D',prsLHS];
        end

    end

    methods (Access = private)
    
        function init(obj, cParams)
            obj.dt            = cParams.dt;
            obj.mesh          = cParams.mesh;
            obj.pressureField = cParams.pressureField;
            obj.velocityField = cParams.velocityField;
        end

        function velLHS = computeVelocityLHS(obj)
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix();
            A = obj.Kelem + obj.Melem;
            s.dim          = obj.velocityField.dim;
            s.globalConnec = obj.velocityField.connec;
            s.nnodeEl      = obj.velocityField.dim.nnodeElem;
            assembler = Assembler(s);
            lhs = assembler.assemble(A);
            velLHS = obj.symGradient(lhs);
        end

        function D = computeDmatrix(obj)
            obj.computeDelem();
            s.dim           = [];
            s.nnodeEl       = [];
            s.globalConnec  = [];
            assembler = Assembler(s);
            D = assembler.assembleFields(obj.Delem, ...
                          obj.velocityField, obj.pressureField);
            obj.D = D;
        end

        function BB = computePressureLHS(obj)
            sz = size(obj.D, 2);
            BB = sparse(sz,sz);
        end

        function A = symGradient(obj, B)
            A = 1/2 * (B+B');
        end

        function computeStiffnessMatrix(obj)
            vel = obj.velocityField;
            ndofs = vel.dim.ndofsElem;
            nelem = obj.mesh.nelem;

%             material = obj.material;
%             Cmat = material.mu;
            
            geom = vel.geometry;
            shape = geom.dNdx;
            ngaus = size(shape,4);
            dvolu = geom.dvolu;
            lhs = zeros(ndofs, ndofs, nelem);
            for igaus = 1:ngaus
                dNdx = shape(:,:,:,igaus);
                dV(1,1,:) = dvolu(:,igaus);
                Bmat = obj.computeB(dNdx);
                Bt   = permute(Bmat,[2 1 3]);
%                 BtC  = pagemtimes(Bt,Cmat);
%                 BtCB = pagemtimes(BtC, Bmat);
                BtB = pagemtimes(Bt,Bmat);
                lhs = lhs + bsxfun(@times, BtB, dV);
            end
            obj.Kelem = lhs;
        end
        
        function computeMassMatrix(obj)
            vel = obj.velocityField;
            nunkn = vel.dim.ndimf;
            nnode = vel.dim.nnodeElem;
            ndofs = vel.dim.ndofsElem;
            nelem = obj.mesh.nelem;
            dtime = obj.dt;
            shpeV = vel.interpolation.shape;
            dvolV = vel.geometry.dvolu;
            ngaus = size(dvolV,2);
            M = zeros(ndofs, ndofs, nelem);
            
            for igauss = 1 :ngaus
                for inode= 1:nnode
                    for jnode= 1:nnode
                        for iunkn= 1:nunkn
                            for junkn= 1:nunkn
                                idof = nunkn*(inode-1)+iunkn;
                                jdof = nunkn*(jnode-1)+junkn;
                                dvol = dvolV(:,igauss);
                                Ni = shpeV(inode,igauss,:);
                                Nj = shpeV(jnode,igauss,:);
                                v = squeeze(Ni.*Nj);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:)/dtime.*dvol;
                            end
                        end
                    end
                end
            end
            obj.Melem = M;
        end
        
        function computeDelem(obj)
            nelem = obj.mesh.nelem;
            vel = obj.velocityField;
            prs = obj.pressureField;
            nunknV = vel.dim.ndimf;
            nnodeV = vel.dim.nnodeElem;
            nnodeP = prs.dim.nnodeElem;
            
            dNdxV = vel.geometry.dNdx;
            dvolV = vel.geometry.dvolu;
            shpeP = prs.interpolation.shape; %nope, should be Quadratic
            ngaus = size(dNdxV,4);

            D = zeros(nunknV*nnodeV,nnodeP,nelem);
            for igauss=1:ngaus
                for inode_var = 1:nnodeP
                    for inode_test = 1:nnodeV
                        for idime = 1:vel.interpolation.ndime
                            dof_test = inode_test*nunknV - nunknV + idime;
                            v = squeeze(dNdxV(idime,inode_test,:,igauss));
                            D(dof_test,inode_var,:)= squeeze(D(dof_test,inode_var,:)) - v(:).*shpeP(inode_var,igauss)...
                                .*dvolV(:,igauss);
                        end
                    end
                end
            end
            obj.Delem = D;
        end
        
        function B = computeB(obj,dNdx)
            vel = obj.velocityField;
            nunkn = vel.dim.ndimf;
            nnode = vel.dim.nnodeElem;
            nelem = obj.mesh.nelem;
            B = zeros(4,nnode*nunkn,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)  = dNdx(1,i,:);
                B(2,j+1,:)= dNdx(1,i,:);
                B(3,j,:)  = dNdx(2,i,:);
                B(4,j+1,:)= dNdx(2,i,:);
            end
        end

    end

end