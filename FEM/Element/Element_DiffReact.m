classdef Element_DiffReact < Element
    %Element_DiffReact Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
        mesh
    end
    
    methods %(Access = ?Physical_Problem)
        function obj = Element_DiffReact(mesh)
            obj.mesh = mesh;
        end
        
        function [r,dr] = computeResidual(obj,uL)
            % *************************************************************
            % Compute
            % - residual: r = Ku - F
            % - residual derivative: dr = K
            % *************************************************************
            [K] = obj.computeStiffnessMatrix;
            [M] = obj.computeMassMatrix(2);
            
            Fext = obj.computeExternalForces();
            R = obj.compute_imposed_displacemet_force(K);
            fext = Fext + R;
            
            fint = K(obj.dof.vF,obj.dof.vF)*uL;
            r = fint - fext;
            dr = epsilon^2*K(obj.dof.vF,obj.dof.vF) + M(obj.dof.vF,obj.dof.vF);
        end
        
        
        function [K] = computeStiffnessMatrix(obj)
            [K] = compute_elem_StiffnessMatrix(obj);
            [K] = obj.AssembleMatrix(K);
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
        
        function [M] = computeMassMatrix(obj,job)
            switch obj.geometry.type
                case 'TRIANGLE'
                    obj.mesh.geometryType = 'Triangle_Linear_Mass';
                case 'QUADRILATERAL'
                    obj.mesh.geometryType = 'Quad_Mass';
            end
            geom = Geometry(obj.mesh);
            dirichlet_data = obj.mesh.connec';
            emat = zeros(geom.nnode,geom.nnode,obj.mesh.nelem);
            
            for igaus=1:geom.ngaus
                for inode=1:geom.nnode
                    for jnode=1:geom.nnode
                        emat(inode,jnode,:)=squeeze(emat(inode,jnode,:)) + geom.weigp(igaus)*geom.shape(inode,igaus)*geom.shape(jnode,igaus)*geom.djacb(:,igaus);
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
                            elumped(inode,:)=elumped(inode,:)+squeeze(emat(inode,jnode,:))';
                        end
                    end
                elseif (nproc==2)
                    for inode=1:nnode
                        for jnode=1:nnode
                            elumped(inode,:)=elumped(inode,:)+squeeze(emat(inode,jnode,:))';
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
                        vmass = squeeze(emat(k,l,:));
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


