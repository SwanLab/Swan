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
    end
    
    methods %(Access = ?Physical_Problem)
        function obj = Element_DiffReact(mesh)
            obj.mesh = mesh;
        end
        function obj = setEpsilon(obj,epsilon)
            obj.epsilon = epsilon;
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
            [M] = obj.AssembleMatrix(M,1,1); % !! UNCOMMENT WHEN INTEGRATION IMPLEMENTED !!
        end
        
        function [K] = compute_elem_StiffnessMatrix(obj)
            
            % Stiffness matrix
            Ke = zeros(obj.nnode,obj.nnode,obj.nelem);
            
            for igaus = 1 :obj.geometry.quadrature.ngaus
                % Strain-displacement matrix
                Bmat = obj.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                
                % Compute Ke
                for iv = 1:obj.nnode
                    for jv = 1:obj.nnode
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
        
%         function compute_elem_MassMatrix(obj)
%             
%         end
        
        %% !! PENDING OF INTEGRATION / ELEMENT DEGREE FOR IMPLEMENTING LIKE STIFFNESS MATRIX !! 
        function [M] = compute_elem_MassMatrix(obj,job)            
            obj.geometry.computeGeometry('QUADRATIC');
%             dirichlet_data = obj.mesh.connec';
            Me = zeros(obj.geometry.interpolation.isoparametric.nnode,obj.geometry.interpolation.isoparametric.nnode,obj.mesh.nelem);
            
            for igaus=1:obj.geometry.quadrature.ngaus
                for inode=1:obj.geometry.interpolation.isoparametric.nnode
                    for jnode=1:obj.geometry.interpolation.isoparametric.nnode
                        Me(inode,jnode,:)=squeeze(Me(inode,jnode,:)) + 0.5*obj.geometry.quadrature.weigp(igaus)*obj.geometry.shape(inode,igaus)*obj.geometry.shape(jnode,igaus)*obj.geometry.djacob(:,igaus);
                    end
                end
            end
            
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
            obj.geometry.computeGeometry('LINEAR');
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
            FextSuperficial = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
    end
    
end


