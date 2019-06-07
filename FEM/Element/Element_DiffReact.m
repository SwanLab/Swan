classdef Element_DiffReact < Element
    %Element_DiffReact Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
        mesh
        K
        M
        Mr
        epsilon
        interpolation_u
    end
    
    properties (Access = private)
        nstre
        addRobinTerm
    end
    
    methods %(Access = ?Physical_Problem)
        function obj = Element_DiffReact(mesh,geometry,material,dof,scale,addRobinTerm)
            obj.mesh = mesh;
            obj.addRobinTerm = addRobinTerm;
            obj.initElement(geometry,material,dof,scale);
            obj.nstre = 2;
            obj.nfields = 1;
            obj.interpolation_u=Interpolation.create(mesh,'LINEAR');
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix(2);
            obj.computeBoundaryMassMatrix();
        end
        
        function obj = setEpsilon(obj,epsilon)
            obj.epsilon = epsilon;
        end
        
        function LHS = computeLHS(obj)
            if obj.addRobinTerm
                LHS = obj.epsilon^2*obj.K + obj.M + 1/obj.epsilon*obj.Mr;              
            else
                LHS = obj.epsilon^2*obj.K + obj.M;
                LHS = obj.bcApplier.full_matrix_2_reduced_matrix(LHS);
            end
        end
        
        function computeStiffnessMatrix(obj)
            Ke = compute_elem_StiffnessMatrix(obj);
            Kg = obj.AssembleMatrix(Ke,1,1); % !!
            obj.K = Kg;
        end
        
        function computeMassMatrix(obj,job)
            Me = compute_elem_MassMatrix(obj,job);
            Mg = obj.AssembleMatrix(Me,1,1); % !!
            obj.M = Mg;
        end
        
        function computeBoundaryMassMatrix(obj)
            if obj.addRobinTerm
                meshB = obj.mesh;
                int = Interpolation.create(meshB,'LINEAR');
                meshType = 'BOUNDARY';
                meshIncludeBoxContour = true;
                cParams = SettingsMeshUnfitted(meshType,meshB,int,meshIncludeBoxContour);
                levelSet = -ones(size(obj.mesh.coord,1),1);
                uMesh = Mesh_Unfitted.create2(cParams);
                uMesh.computeMesh(levelSet);
                integrator = Integrator.create(uMesh);
                obj.Mr = integrator.integrateLHS(uMesh);
            end
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
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
            % !!!!!!!!!!!!!!!!!!!!
            
            M = Me;
            
        end
    end
    
    methods(Access = protected) % Only the child sees the function
        function FextSuperficial = computeSuperficialFext(obj)
            FextSuperficial = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj)
            FextVolumetric = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
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
end


