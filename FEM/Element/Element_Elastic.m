classdef Element_Elastic < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fext
        interpolation_u
        K
        K_generator
        DeltaC
    end
    
    methods (Static) %(Access = {?Physical_Problem, ?Element_Elastic_Micro, ?obj})
        function obj = create(mesh,geometry,material,dof)
            switch mesh.scale
                case 'MICRO'
                    switch mesh.pdim
                        case '2D'
                            obj = Element_Elastic_2D_Micro(mesh,geometry,material,dof);
                        case '3D'
                            obj = Element_Elastic_3D_Micro(mesh,geometry,material,dof);
                    end
                case 'MACRO'
                    switch mesh.pdim
                        case '2D'
                            obj = Element_Elastic_2D(mesh,geometry,material,dof);
                        case '3D'
                            obj = Element_Elastic_3D(mesh,geometry,material,dof);
                    end
            end
        end
    end
    
    methods %(Access = {?Physical_Problem, ?Element_Elastic_Micro, ?Element})
        function obj = Element_Elastic(mesh,geometry,material,dof,nstre)
            obj@Element(geometry,material,dof)
            obj.nstre = nstre;
            obj.nfields=1;
            obj.interpolation_u = Interpolation.create(mesh,'LINEAR');
            
            obj.initialize_dvolum()
            Bmat   = obj.computeBmat();
            dvolum = obj.geometry.dvolu;
            ngaus  = obj.quadrature.ngaus;
            
            dim                = DimensionVariables();
            dim.nnode          = obj.nnode;
            dim.nunkn          = obj.dof.nunkn;
            dim.nstre          = obj.nstre;
            dim.ndof           = obj.dof.ndof;
            dim.nelem          = obj.nelem;
            dim.ndofPerElement = dim.nnode*dim.nunkn;
            dim.ngaus          = ngaus;
            dim.nentries       = dim.nelem*(dim.ndofPerElement)^2;
            
            connec = obj.geometry.interpolation.T;
            
%            obj.DeltaC = ContitutiveTensorIncrement();
            
            obj.K_generator = StiffnessMatrixGenerator(connec,Bmat,dvolum,dim);
        end
        
        function Bmat = computeBmat(obj)
            ngaus = obj.quadrature.ngaus;
            nnode = obj.interpolation_u.nnode;
            nunkn = obj.dof.nunkn;
            nstre = obj.nstre;
            nelem = obj.nelem;
            
            Bmat = zeros(ngaus,nstre,nnode*nunkn,nelem);
            for igaus = 1:ngaus
                Bmat(igaus,:,:,:) = obj.computeB(igaus);
            end
        end
        
        
        function initialize_dvolum(obj)
            obj.computeQuadrature()
            obj.computeInterpolation()
            obj.computeGeometry() 
        end
        
        function computeQuadrature(obj)
            obj.quadrature.computeQuadrature('LINEAR');
        end
        
        function computeInterpolation(obj)
            obj.interpolation_u.computeShapeDeriv(obj.quadrature.posgp)
        end
        
        function computeGeometry(obj)
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
        end
        
        function Kred = computeLHS(obj)
            obj.K = obj.computeStiffnessMatrix;
            Kred = obj.full_matrix_2_reduced_matrix(obj.K);
        end
        
        function fext_red = computeRHS(obj)
            Fext = obj.computeExternalForces;
            R = obj.compute_imposed_displacement_force(obj.K);
            obj.fext = Fext + R;
            fext_red = obj.full_vector_2_reduced_vector(obj.fext);
        end
        
        function [K] = computeStiffnessMatrix(obj)
            [K] = obj.computeStiffnessMatrixSYM;            
        end
        
        function [idx1,idx2,nunkn1,nunkn2,nnode1,nnode2,col,row] = get_assemble_parameters(obj,~,~)
            idx1 = cell2mat(obj.dof.in_elem);
            idx2 = idx1;
            nunkn1 = obj.dof.nunkn;
            nnode1 = obj.interpolation_u.nnode;
            nunkn2 = nunkn1;
            nnode2 = nnode1;
            col = obj.dof.ndof;
            row = col;
        end
        
        function K = compute_elem_StiffnessMatrix(obj)
            obj.quadrature.computeQuadrature('LINEAR');
            obj.interpolation_u.computeShapeDeriv(obj.quadrature.posgp)
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
            % Stiffness matrix
            Ke = zeros(obj.dof.nunkn*obj.nnode,obj.dof.nunkn*obj.nnode,obj.nelem);
            
            % Elastic matrix
            Cmat = obj.material.C;
            for igaus = 1 :obj.quadrature.ngaus
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
        function K = computeStiffnessMatrixSYM(obj)
           %  obj.DeltaC.obtainChangedElements(obj.material.C)
             obj.K_generator.generate(obj.material.C);
             K = obj.K_generator.getStiffMatrix();
             
             
%              ndofPerElement = obj.nnode*obj.dof.nunkn;
%              Bfull = zeros(obj.quadrature.ngaus,obj.nstre,ndofPerElement,obj.nelem);
%              Bfull2 = zeros(obj.quadrature.ngaus*obj.nelem*obj.nstre,ndofPerElement);
%              for igaus = 1:obj.quadrature.ngaus
%                  unitaryIndex = false(obj.quadrature.ngaus*obj.nstre,1);
%                  pos = obj.nstre*(igaus-1) + 1 : obj.nstre*(igaus) ;
%                  unitaryIndex(pos) = true;
%                  
%                 Index = repmat(unitaryIndex,obj.nelem,1); 
%                 Bfull(igaus,:,:,:) = obj.computeB(igaus);
%                 Bshif = reshape(permute(obj.computeB(igaus),[1 3 2]),obj.nelem*obj.nstre,ndofPerElement);
%                 Bfull2(Index,:) = Bshif;
%              end
%              
%             dim                = DimensionVariables();
%             dim.nnode          = obj.nnode;
%             dim.nunkn          = obj.dof.nunkn;
%             dim.nstre          = obj.nstre;
%             dim.ndof           = obj.dof.ndof;
%             dim.nelem          = obj.nelem;
%             dim.ndofPerElement = dim.nnode*dim.nunkn;
%             dim.ngaus          = obj.quadrature.ngaus;
%             dim.nentries       = dim.nelem*(dim.ndofPerElement)^2;
%              
%              conec = obj.geometry.interpolation.T;
%              K2 = KGeneratorWithfullStoredB(Bfull2,dim,conec,obj.material.C,obj.geometry.dvolu);
%              K = K2.K; 
        end
        
        
    end
    
    methods(Access = protected) 
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
            variables.stress = obj.computeStress(variables.strain,obj.material.C,obj.quadrature.ngaus,obj.nstre);
        end
        
        function u = compute_displacements(obj,usol)
            u = obj.reduced_vector_2_full_vector(usol);
        end
        
        function strain = computeStrain(obj,d_u,idx)
            strain = zeros(obj.nstre,obj.nelem,obj.quadrature.ngaus);
            for igaus = 1:obj.quadrature.ngaus
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
