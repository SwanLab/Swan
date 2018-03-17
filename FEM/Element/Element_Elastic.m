classdef Element_Elastic < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fext
        K
    end
    
    methods (Static) %(Access = {?Physical_Problem, ?Element_Elastic_Micro, ?obj})
        function obj = create(mesh,geometry,material,dof)
            switch mesh.scale
                case 'MICRO'
                    obj = Element_Elastic_2D_Micro(geometry,material,dof);
                    obj.nstre = 3;
                case 'MACRO'
                    switch mesh.pdim
                        case '2D'
                            obj = Element_Elastic_2D(geometry,material,dof);
                            obj.nstre = 3;
                        case '3D'
                            obj = Element_Elastic_3D(geometry,material,dof);
                            obj.nstre = 6;
                    end
            end
        end
    end
    
    methods %(Access = ?Elastic_Problem)
        function obj = Element_Elastic(geometry,material,dof)
            obj = obj@Element(geometry,material,dof);
        end
        
        function r = computeResidual(obj,x,Kred)
            
            %             [K] = obj.computeStiffnessMatrix();
            
            Fext = obj.computeExternalForces();
            R = obj.compute_imposed_displacemet_force(obj.K);
            obj.fext = Fext + R;
            
            %             Kred = obj.full_matrix_2_reduced_matrix(K);
            fext_red = obj.full_vector_2_reduced_vector(obj.fext);
            
            fint_red = Kred*x;
            
            r = fint_red - (fext_red);
            %             dr = Kred;
        end
        
        function [K] = computeStiffnessMatrix(obj)
            K = compute_elem_StiffnessMatrix(obj);
            [K] = obj.AssembleMatrix(K,1,1);
        end
        
        function dr = computedr(obj)
            obj.K = obj.computeStiffnessMatrix;
            Kred = obj.full_matrix_2_reduced_matrix(obj.K);
            dr = Kred;
        end
        
        function K = compute_elem_StiffnessMatrix(obj)
            % Stiffness matrix
            Ke = zeros(obj.dof.nunkn*obj.nnode,obj.dof.nunkn*obj.nnode,obj.nelem);
            
            % Elastic matrix
            Cmat = obj.material.C;
            for igaus = 1 :obj.geometry.quadrature.ngaus
                % Strain-displacement matrix
                Bmat = obj.computeB(igaus);
                
                for iv = 1:obj.nnode*obj.dof.nunkn
                    for jv = 1:obj.nnode*obj.dof.nunkn
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
        function FextSuperficial = computeSuperficialFext(obj)
            FextSuperficial = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj)
            FextVolumetric = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
        
        function variables = computeDispStressStrain(obj,uL)
            variables.d_u = obj.compute_displacements(uL);
            variables.fext = obj.fext;
            variables.strain = obj.computeStrain(variables.d_u,obj.dof.in_elem{1});
            variables.stress = obj.computeStress(variables.strain,obj.material.C,obj.geometry.quadrature.ngaus,obj.nstre);
        end
        
        function u = compute_displacements(obj,usol)
            u = obj.reduced_vector_2_full_vector(usol);
        end
        
        function strain = computeStrain(obj,d_u,idx)
            strain = zeros(obj.nstre,obj.nelem,obj.geometry.quadrature.ngaus);
            for igaus = 1:obj.geometry.quadrature.ngaus
                Bmat = obj.computeB(igaus);
                %Bmat = Bmat{1,1};
                for istre=1:obj.nstre
                    for inode=1:obj.nnode
                        for idime=1:obj.dof.nunkn
                            ievab = obj.dof.nunkn*(inode-1)+idime;
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
    end
    
    methods(Static, Access = protected)
        % Only used in Element_Elastic_2D
        function strain = computeEz(strain,nstre,nelem,material)
            mu = material.mu;
            kappa = material.kappa;
            epoiss = (kappa(1,1) - mu(1,1))./(kappa(1,1) + mu(1,1));
            epoiss = full(ones(1,nelem)*epoiss);
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
