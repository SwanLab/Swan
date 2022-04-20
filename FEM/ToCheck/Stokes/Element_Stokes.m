classdef Element_Stokes < Element
    %Element_Stokes Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LHS_elem
        LHS
        RHS
        interpolation_v
        interpolation_p

        dim
    end

    properties(Access = private)
        D_elem
        M_elem
        K_elem
        D
        dt
        mesh
    end
    
    methods
        function obj = Element_Stokes(geometry,mesh,material,dof,problemData,interp, dim)
            obj.initElement(geometry,mesh,material,dof,problemData.scale,interp);
            obj.dim = dim;
            obj.mesh = mesh;
            %obj.nstre=0;
            obj.nfields=2;
            obj.interpolation_v= interp{1};
            obj.interpolation_p= interp{2};
        end
        
        function [r,dr] = computeResidual(obj,x,dr,x_n)
            %             K = compute_LHS(obj);
            if (nargin ==3)
                % Steady
                Mred_x_n = zeros(length(obj.dof.free{1}),1);
            else
                % Transient
                M = obj.AssembleMatrix(obj.M_elem,1,1);
                M = obj.symGradient(M);
                Mred = M(obj.dof.free{1},obj.dof.free{1});
                Mred_x_n = Mred*x_n;
            end
            
            Fext = compute_RHS(obj);
            
            
            R = obj.compute_imposed_displacement_force(obj.LHS);
            Fext = Fext + R ;
            
            
            Fext_red = obj.bcApplier.fullToReducedVector(Fext);
            Fext_red(1:length(obj.dof.free{1}),1) = Fext_red(1:length(obj.dof.free{1}),1) + Mred_x_n;
            
            fint_red = dr*x;
            
            r = fint_red - (Fext_red);
            %             dr = Kred;
            
        end
        
        function dr = computedr(obj,dt)
            if nargin < 2
                dt=inf;
            end
            obj.LHS = compute_LHS(obj,dt);
            LHSred = obj.bcApplier.fullToReducedMatrix(obj.LHS);
            dr = LHSred;
        end
        
        function LHS = compute_LHS(obj,dt)
            obj.dt = dt;
            AA = obj.computeAAmatrix();
            D = obj.computeDmatrix();
            BB = obj.computeBBmatrix();
            LHS = [AA, D; D',BB];
        end
        
        function RHS = compute_RHS(obj)
            Fext = obj.computeVolumetricFext(obj.nelem,obj.geometry,obj.dof);
            g = obj.compute_velocity_divergence;
            RHS_elem{1,1} = Fext;
            RHS_elem{2,1} = g;
            RHS = AssembleVector(obj,RHS_elem);
        end
        
        function Fext = compute_vol_force_on_nodes(obj,geometry,idx,nnode,nunkn)
            %             for i = 1:length(bc.iN)
            %                 for j = 1:nnode*nunkn
            %                     ind = find(idx(j,:) == bc.iN(i));
            %                     if ~isempty(ind)
            %                         f(j,:,ind) = bc.force(i,3);
            %                     end
            %                     %                     clear ind
            %                     ind = [];
            %                 end
            %             end
            
            
            %              for igaus=1:geometry.ngaus
            %                 for inode=1:nnode
            %                     for jnode=1:nnode
            %                         for iunkn=1:nunkn
            %                             elemental_dof = jnode*nunkn-nunkn+iunkn; %% dof per guardar el valor de la integral
            %
            %                                 v= squeeze(geometry.shape(inode,igaus).*geometry.shape(jnode,igaus).*f(elemental_dof,1,:));
            %                                 Fext(elemental_dof,1,:)= squeeze(Fext(elemental_dof,1,:)) + v(:).*geometry.dvolu(:,igaus);
            %
            %                         end
            %                     end
            %                 end
            %             end
        end
        
        function Fext = compute_vol_force_on_gauss_points(obj,geometry,nnode,nunkn,f)
            Fext = zeros(nnode*nunkn,1,obj.nelem);
            for igaus=1:obj.quadrature.ngaus
                for inode=1:nnode
                    for iunkn=1:nunkn
                        elemental_dof = inode*nunkn-nunkn+iunkn; %% dof per guardar el valor de la integral
                        shape = obj.interpolation_v.shape(inode,igaus);
                        fvalue = f(iunkn,igaus,:);
                        v= squeeze(shape.*fvalue);
                        Fext(elemental_dof,1,:)= squeeze(Fext(elemental_dof,1,:)) + v(:).*geometry.dvolu(:,igaus);
                        
                    end
                end
            end
        end
        
        function g = compute_velocity_divergence(obj)
            dimP = obj.dim{1};
            nunkn = dimP.ndimField;
            g = zeros(obj.interp{2}.nnode*nunkn,1,obj.nelem);
        end
        
        function variable = computeVars(obj,x_free)
            x = obj.bcApplier.reducedToFullVector(x_free);
            variable.u = x(1:obj.dof.ndof(1),:);
            variable.p = x(obj.dof.ndof(1)+1:end,:);
        end
    end
    
    methods (Access = protected)
        
        function Fext = computePuntualRHS(obj,nunkn,nelem,nnode,bc,idx)
            Fext = zeros(nnode*nunkn,1,nelem);
            for i = 1:length(bc.iN)
                for j = 1:nelem
                    ind = find(idx(:,j) == bc.iN(i));
                    if ~isempty(ind)
                        Fext(ind,:,j) = bc.neunodes(i,3);
                    end
                    % clear ind
                    ind = [];
                end
            end
        end
        
        function Fext = computeSuperficialFext(obj,nunkn,nelem,nnode,bc,idx) %To be donne
            % Fext = zeros(nnode*nunkn,1,nelem);
            Fext = 0;
        end
        
        function Fext = computeVolumetricFext(obj,nelem,geometry,dof)
            idx = obj.dof.in_elem{1};
            geometry = geometry(1);
            nnode = obj.interpolation_v.nnode;
            dimV = obj.dim{1};
            nunkn = dimV.ndimField;
            %             f = zeros(nnode*nunkn,1,nelem);
            
            %             obj.RHS = zeros(nnode*nunkn,1,nelem);
            
            if  ~isempty(dof.neumann_values)
                if ~ismatrix(dof.neumann_values)
                    f=dof.neumann_values;
                    Fext = obj.compute_vol_force_on_gauss_points(geometry,nnode,nunkn,f);
                else
                    Fext = obj.compute_vol_force_on_nodes(geometry,idx,nnode,nunkn);
                end
            else
                Fext = zeros(nnode*nunkn,1,nelem);
            end
        end
    end

    methods (Access = private)

        function A = symGradient(obj, B)
            A = 1/2 * (B+B');
        end

        function Aassembled = computeAAmatrix(obj)
            obj.compute_M();
            obj.compute_K();
            A = obj.K_elem + obj.M_elem;
            Aass = obj.AssembleMatrix(A ,1, 1);
            Aassembled = obj.symGradient(Aass);
        end

        function Dassembled = computeDmatrix(obj)
            obj.compute_D();
            Dassembled = obj.AssembleMatrix(obj.D_elem, 1, 2);
            obj.D = Dassembled;
        end

        function BB = computeBBmatrix(obj)
            sz = size(obj.D, 2);
            BB = sparse(sz,sz);
        end
        
        function M = compute_M(obj)
            dimV = obj.dim{1};
            nunkn = dimV.ndimField;
            nnode = obj.interpolation_v.nnode;
            ndofs = nunkn*nnode;
            nelem = obj.nelem;
            dtime = obj.dt;
            M = zeros(ndofs, ndofs, nelem);
            
            for igauss = 1 :obj.quadrature.ngaus
                for inode= 1:nnode
                    for jnode= 1:nnode
                        for iunkn= 1:nunkn
                            for junkn= 1:nunkn
                                idof = nunkn*(inode-1)+iunkn;
                                jdof = nunkn*(jnode-1)+junkn;
                                dvol = obj.geometry(1).dvolu(:,igauss);
                                Ni = obj.interpolation_v.shape(inode,igauss,:);
                                Nj = obj.interpolation_v.shape(jnode,igauss,:);
                                v = squeeze(Ni.*Nj);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:)/dtime.*dvol;
                            end
                        end
                    end
                end
            end
            obj.M_elem = M;
%             obj.computeMassMatrix();
        end

        function computeMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS'; % INTERPOLATIONTYPE
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim{1};
            LHS = LHSintegrator.create(s);
            Mass = LHS.compute();
        end

        function computeStiffnessMatrix(obj)
            s.type = 'ElasticStiffnessMatrix';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim{1};
            s.material     = obj.material;
%             s.material.C = obj.material.mu;
            LHS = LHSintegrator.create(s);
            K   = LHS.compute();
            Kred = obj.boundaryConditions.fullToReducedMatrix(K);
            obj.stiffnessMatrix    = K;
            obj.stiffnessMatrixRed = Kred;
        end
        
        function K = compute_K(obj)
            dimV = obj.dim{1};
            nunkn = dimV.ndimField;
            nnode = obj.interpolation_v.nnode;
            ndofs = nunkn*nnode;
            nelem = obj.nelem;
            material = obj.material;
            K = zeros(ndofs, ndofs, nelem);
            
            Cmat = material.mu;
            obj.quadrature.computeQuadrature('QUADRATIC');
            obj.geometry(1).computeGeometry(obj.quadrature,obj.interpolation_v);
            for igauss = 1 :obj.quadrature.ngaus
                dNdx = obj.geometry(1).dNdx(:,:,:,igauss);
                Bmat = obj.computeB(nunkn, nelem, nnode, dNdx);
                for iv=1:ndofs
                    for jv=1:ndofs
                        for istre=1:nunkn*obj.interpolation_v.ndime
                            for jstre=1:nunkn*obj.interpolation_v.ndime
                                v = squeeze(Bmat(istre,iv,:).*Cmat(istre,jstre,:).*Bmat(jstre,jv,:));
                                K(iv,jv,:)= squeeze(K(iv,jv,:)) + v(:).*obj.geometry(1).dvolu(:,igauss);
                                %                                 obj.LHS(iv,jv,:) = squeeze(obj.LHS(iv,jv,:)) + v(:).*geometry.dvolu(:,igauss);
                            end
                        end
                    end
                end
                
            end
            obj.K_elem = K;
%             obj.computeStiffnessMatrix();
        end
        
        function D = compute_D(obj)
            dimV = obj.dim{1};
            nelem = obj.nelem;
            nunknU = dimV.ndimField;
            nnodeV = obj.interpolation_v.nnode;
            
            D = zeros(nunknU*obj.interpolation_v.nnode,obj.interpolation_p.nnode,nelem);
            obj.quadrature.computeQuadrature('QUADRATIC');
            obj.geometry(2).computeGeometry(obj.quadrature,obj.interpolation_p);
            for igauss=1:obj.quadrature.ngaus
                for inode_var = 1:obj.interpolation_p.nnode
                    for inode_test = 1:nnodeV
                        for idime = 1:obj.interpolation_v.ndime
                            dof_test = inode_test*nunknU - nunknU + idime;
                            v= squeeze (obj.geometry(1).dNdx(idime,inode_test,:,igauss));
                            D(dof_test,inode_var,:)= squeeze(D(dof_test,inode_var,:)) - v(:).*obj.interpolation_p.shape(inode_var,igauss)...
                                .*obj.geometry(1).dvolu(:,igauss);
                        end
                    end
                end
            end
            obj.D_elem = D;
        end
        
        function B = computeB(obj,nunkn,nelem,nnode,dNdx)
            B = zeros(2,nnode*nunkn,nelem);
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