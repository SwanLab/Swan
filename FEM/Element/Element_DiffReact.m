classdef Element_DiffReact < Element
    %Element_DiffReact Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
        mesh
        K
        M
        epsilon
        interpolation_u
    end
    
    methods %(Access = ?Physical_Problem)
        function obj = setEpsilon(obj,epsilon)
            obj.epsilon = epsilon;
        end
        function obj=Element_DiffReact(mesh,geometry,material,dof)
            obj@Element(geometry,material,dof);
            obj.nstre=2;
            obj.nfields=1;
            obj.interpolation_u=Interpolation.create(mesh,'LINEAR');
        end
        function [r,dr] = computeResidual(obj,uL)
            % *************************************************************
            % Compute
            % - residual: r = Ku - F
            % - residual derivative: dr = K
            % *************************************************************
            [K] = obj.computeStiffnessMatrix;
            [M] = obj.computeMassMatrix(2);
            obj.K = K; obj.M = M;
            
            Fext = obj.computeExternalForces();
            R = obj.compute_imposed_displacemet_force(K);
            fext = Fext + R;
            
            Kred = K(obj.dof.free,obj.dof.free);
            Mred = M(obj.dof.free,obj.dof.free);
            
            fint = Kred*uL;
            r = fint - fext;
            dr = obj.epsilon^2*Kred + Mred;
        end
        
        
        function [K] = computeStiffnessMatrix(obj)
            [K] = compute_elem_StiffnessMatrix(obj);
            [K] = obj.AssembleMatrix(K,1,1);
        end
        
        function [M] = computeMassMatrix(obj,job)
            [M] = compute_elem_MassMatrix(obj,job);
            [M] = obj.AssembleMatrix(M,1,1);
        end
        
        function [K] = compute_elem_StiffnessMatrix(obj)
            obj.quadrature.computeQuadrature('LINEAR');
            obj.interpolation_u.computeShapeDeriv(obj.quadrature.posgp)
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
            % Stiffness matrix
            Ke = zeros(obj.dof.nunkn*obj.nnode,obj.dof.nunkn*obj.nnode,obj.nelem);
            
            for igaus = 1 :obj.quadrature.ngaus
                % Strain-displacement matrix
                Bmat = obj.computeB(obj.dof.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                
                % Compute Ke
                
                for iv = 1:obj.nnode*obj.dof.nunkn
                    for jv = 1:obj.nnode*obj.dof.nunkn
                        for istre = 1:obj.nstre
                            % for jstre=1:nstre
                            v = squeeze(Bmat(istre,iv,:).*Bmat(istre,jv,:));
                            Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*obj.geometry.dvolu(:,igaus);
                            %end
                        end
                    end
                end
            end
            
            K = Ke;
        end

        function [M] = compute_elem_MassMatrix(obj,job)            
            obj.quadrature.computeQuadrature('QUADRATICMASS');
            obj.interpolation_u.computeShapeDeriv(obj.quadrature.posgp)
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
            Me = zeros(obj.interpolation_u.nnode,obj.interpolation_u.nnode,obj.nelem);
            
            for igaus=1:obj.quadrature.ngaus
                for inode=1:obj.interpolation_u.nnode
                    for jnode=1:obj.interpolation_u.nnode
                        Me(inode,jnode,:)=squeeze(Me(inode,jnode,:)) + obj.quadrature.weigp(igaus)*obj.interpolation_u.shape(inode,igaus)...
                            *obj.interpolation_u.shape(jnode,igaus)*obj.geometry.djacob(:,igaus);
                    end
                end
            end
            %% PENDING TO BE REMOVED AS SOON AS PHYSPROBLEM FAMILIY IS IMPLEMENTED
            obj.quadrature.computeQuadrature('LINEAR');
            obj.interpolation_u.computeShapeDeriv(obj.quadrature.posgp)
            %% !!!!!!!!!!!!!!!!!!!!
            
            M = Me;
            
            %% !! ERROR IN QUADRATURE: NGAUS = 1, WHEN THE
%             
%             if (job==1)
%                 % lumped mass matrix
%                 elumped = zeros(obj.geometry.nnode,obj.mesh.nelem);
%                 M = zeros(obj.geom.nnode,1);
%                 [nproc,coeff] = nprocedure(etype,nnode);
%                 if (nproc==1)
%                     for inode=1:nnode
%                         for jnode=1:nnode
%                             elumped(inode,:)=elumped(inode,:)+squeeze(Me(inode,jnode,:))';
%                         end
%                     end
%                 elseif (nproc==2)
%                     for inode=1:nnode
%                         for jnode=1:nnode
%                             elumped(inode,:)=elumped(inode,:)+squeeze(Me(inode,jnode,:))';
%                         end
%                         elumped(inode,:)=elumped(inode,:)*coeff(inode);
%                     end
%                 end
%                 for inode=1:nnode
%                     M = M + sparse(dirichlet_data(inode,:),1,elumped(inode,:),npnod,1);
%                 end
%             elseif (job==2)
%                 M = sparse(obj.mesh.npnod*obj.geometry.quadrature.ngauss,obj.mesh.npnod*obj.geometry.quadrature.ngaus);
%                 for k=1:obj.geometry.quadrature.ngaus
%                     for l=1:obj.geometry.quadrature.ngaus
%                         vmass = squeeze(Me(k,l,:));
%                         M = M + sparse(dirichlet_data(k,:),dirichlet_data(l,:),vmass,obj.mesh.npnod,obj.mesh.npnod);
%                     end
%                 end
%             end
        end        
    end
    
    methods (Static)
        function [B] = computeB(nunkn,nelem,nnode,cartd)
            B = zeros(2,nnode*nunkn,nelem);
            for inode=1:nnode
                j = nunkn*(inode-1)+1;
                B(1,j,:)=cartd(1,inode,:);
                B(2,j,:)=cartd(2,inode,:);
            end
        end
    end
    
    
    methods (Access = protected)
        function FextSuperficial = computeSuperficialFext(obj,bc)
            FextSuperficial = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
    end
    
end


