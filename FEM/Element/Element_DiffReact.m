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
            [K] = obj.AssembleMatrix(K);
        end
        
        function [M] = computeMassMatrix(obj,job)
            [M] = compute_elem_MassMatrix(obj,job);
%             [M] = obj.AssembleMatrix(M); % !! UNCOMMENT WHEN INTEGRATION IMPLEMENTED !!
        end
        
        function [K] = compute_elem_StiffnessMatrix(obj)
            
            % Stiffness matrix
            Ke = zeros(obj.nunkn*obj.nnode,obj.nunkn*obj.nnode,obj.nelem);
            
            for igaus = 1 :obj.geometry.ngaus
                % Strain-displacement matrix
                Bmat = obj.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                
                % Compute Ke
                if obj.nelem < 1000 %Just to reduce test.m compute time TO BE REMOVED
                    for i = 1:obj.nelem
                        Ke(:,:,i) = Ke(:,:,i)+Bmat(:,:,i)'*Bmat(:,:,i)*obj.geometry.dvolu(i,igaus);
                    end
                else
                    for iv = 1:obj.nnode*obj.nunkn
                        for jv = 1:obj.nnode*obj.nunkn
                            for istre = 1:obj.nstre
                                % for jstre=1:nstre
                                v = squeeze(Bmat(istre,iv,:).*Bmat(istre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*obj.geometry.dvolu(:,igaus);
                                %end
                            end
                        end
                    end
                end
            end
            K = Ke;
        end
        
        %% !! PENDING OF INTEGRATION / ELEMENT DEGREE FOR IMPLEMENTING LIKE STIFFNESS MATRIX !!
        function [M] = compute_elem_MassMatrix(obj,job)
            %             % Stiffness matrix
            %             Me = zeros(obj.nunkn*obj.nnode,obj.nunkn*obj.nnode,obj.nelem);
            %
            %             for igaus = 1 :obj.geometry.ngaus
            %                 % Compute Me
            %                 for iv = 1:obj.nnode*obj.nunkn
            %                     for jv = 1:obj.nnode*obj.nunkn
            %                         for istre = 1:obj.nstre
            %                             % for jstre=1:nstre
            %                             Me(iv,jv,:) = squeeze(Me(iv,jv,:)) + obj.geometry.weigp(igaus)*obj.geometry.shape(iv,igaus)*obj.geometry.shape(jv,igaus)*obj.geometry.djacb(:,igaus);
            %                             %end
            %                         end
            %                     end
            %                 end
            %             end
            %             M = Me;
            %         end
            
            switch obj.geometry.type
                case 'TRIANGLE'
                    obj.mesh.geometryType = 'Triangle_Linear_Mass';
                case 'QUADRILATERAL'
                    obj.mesh.geometryType = 'Quad_Mass';
            end
            
            geom = Geometry(obj.mesh);
            dirichlet_data = obj.mesh.connec';
            Me = zeros(geom.nnode,geom.nnode,obj.mesh.nelem);
            
            for igaus=1:geom.ngaus
                for inode=1:geom.nnode
                    for jnode=1:geom.nnode
                        Me(inode,jnode,:)=squeeze(Me(inode,jnode,:)) + geom.weigp(igaus)*geom.shape(inode,igaus)*geom.shape(jnode,igaus)*geom.djacb(:,igaus);
                    end
                end
            end
            
            if (job==1)
                % lumped mass matrix
                elumped = zeros(obj.geom.nnode,obj.mesh.nelem);
                M = zeros(obj.geom.nnode,1);
                [nproc,coeff] = nprocedure(etype,nnode);
                if (nproc==1)
                    for inode=1:nnode
                        for jnode=1:nnode
                            elumped(inode,:)=elumped(inode,:)+squeeze(Me(inode,jnode,:))';
                        end
                    end
                elseif (nproc==2)
                    for inode=1:nnode
                        for jnode=1:nnode
                            elumped(inode,:)=elumped(inode,:)+squeeze(Me(inode,jnode,:))';
                        end
                        elumped(inode,:)=elumped(inode,:)*coeff(inode);
                    end
                end
                for inode=1:nnode
                    M = M + sparse(dirichlet_data(inode,:),1,elumped(inode,:),npnod,1);
                end
            elseif (job==2)
                M = sparse(obj.mesh.npnod,obj.mesh.npnod);
                for k=1:geom.nnode
                    for l=1:geom.nnode
                        vmass = squeeze(Me(k,l,:));
                        M = M + sparse(dirichlet_data(k,:),dirichlet_data(l,:),vmass,obj.mesh.npnod,obj.mesh.npnod);
                    end
                end
            end
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


