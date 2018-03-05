classdef Element_Elastic < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here

    
    properties
          fext
    end
    
    methods %(Access = {?Physical_Problem, ?Element_Elastic_Micro, ?Element})
        function [r,dr] = computeResidual(obj,x)

            [K] = obj.computeStiffnessMatrix();
            
            Fext = obj.computeExternalForces();            
            R = obj.compute_imposed_displacemet_force(K);
            obj.fext = Fext + R;
            
            Kred = obj.full_matrix_2_reduced_matrix(K,obj.dof);            
            fext_red = obj.full_vector_2_reduced_vector(obj.fext,obj.dof);

            fint_red = Kred*x;

            r = fint_red - (fext_red);
            dr = Kred;
        end
        
        function [K] = computeStiffnessMatrix(obj)
            [K] = compute_elem_StiffnessMatrix(obj);                        
            [K] = obj.AssembleMatrix(K);
        end
        

        

        
        function [K] = compute_elem_StiffnessMatrix(obj)
            
            % Stiffness matrix
            Ke = zeros(obj.nunkn*obj.nnode,obj.nunkn*obj.nnode,obj.nelem);
            
            % Elastic matrix
            Cmat = obj.material.C;
            for igaus = 1 :obj.geometry.ngaus
                % Strain-displacement matrix
                Bmat = obj.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                
                for iv = 1:obj.nnode*obj.nunkn
                    for jv = 1:obj.nnode*obj.nunkn
                        for istre = 1:obj.nstre
                            for jstre = 1:obj.nstre
                                v = squeeze(Bmat(istre,iv,:).*Cmat(istre,jstre,:).*Bmat(jstre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*obj.geometry.dvolu(:,igaus);
                            end
                        end
                        
                    end
                end
                
            end
            K = Ke;
        end
        
    end
    
    
    methods(Access = protected) % Only the child sees the function
        function FextSuperficial = computeSuperficialFext(obj,bc)
            FextSuperficial = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function variables = computeDispStressStrain(obj,uL)
            variables.d_u = obj.compute_displacements(uL);
            variables.fext = obj.fext;
            variables.strain = obj.computeStrain(variables.d_u,obj.dim,obj.nnode,obj.nelem,obj.geometry.ngaus,obj.dof.in_elem{1});
            variables.stress = obj.computeStress(variables.strain,obj.material.C,obj.geometry.ngaus,obj.nstre);
        end
        
       
        function u = compute_displacements(obj,usol)
            u = obj.reduced_vector_2_full_vector(usol,obj.dof);
        end
        
        
        function strain = computeStrain(obj,d_u,dim,nnode,nelem,ngaus,idx)
            strain = zeros(dim.nstre,nelem,ngaus);
            for igaus = 1:ngaus
                Bmat = obj.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                %Bmat = Bmat{1,1};
                for istre=1:dim.nstre
                    for inode=1:nnode
                        for idime=1:dim.nunkn
                            ievab = dim.nunkn*(inode-1)+idime;
                            strain(istre,:,igaus)=strain(istre,:,igaus)+(squeeze(Bmat(istre,ievab,:)).*d_u(idx(ievab,:)))';
                        end
                    end
                end
            end
        end
        
    end
    
    methods(Static)
        function variables = permuteStressStrain(variables)
            variables.strain = permute(variables.strain, [3 1 2]);
            variables.stress = permute(variables.stress, [3 1 2]);
        end
        
        function Ared = full_matrix_2_reduced_matrix(A,dof)
            Ared = A(dof.free{1},dof.free{1});
        end
        
        function b_red = full_vector_2_reduced_vector(b,dof)
            b_red = b(dof.free{1});
        end
        
        function b = reduced_vector_2_full_vector(bfree,dof)
            b = zeros(dof.ndof,1);
            b(dof.free{1}) = bfree;
            b(dof.dirichlet{1}) = dof.dirichlet_values{1};
        end
        
    end
    
    methods(Static, Access = protected)
        % Only used in Element_Elastic_2D
        function strain = computeEz(strain,nstre,nelem,material)
            mu = material.mu;
            kappa = material.kappa;
            epoiss = (kappa(1,1) - mu(1,1))./(kappa(1,1) + mu(1,1));
            epoiss = ones(1,nelem)*epoiss;
            strain(nstre+1,:,:) = (-epoiss./(1-epoiss)).*(strain(1,:,:)+strain(2,:,:));
        end
        
        % Compute strains (e = B�u)
        
        
        % Compute stresses
        function stres = computeStress(strain,C,ngaus,nstre)
            stres = zeros(size(strain));
            for igaus = 1:ngaus
                for istre=1:nstre
                    for jstre=1:nstre
                        stres(istre,:,igaus) = stres(istre,:,igaus) + squeeze(C(istre,jstre,:))'.*strain(jstre,:,igaus);
                    end
                end
            end
        end
    end
end
