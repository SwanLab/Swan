classdef Element_Elastic < Element

    properties (Access = public)
        fextRed
    end
        
    properties
        fext
        interpolation_u
        K
        K_generator
        DeltaC
        StiffnessMatrix
        dim
        connec
        Bmatrix
    end
    
    properties (Access = protected, Abstract)
        nstre
    end
    
    properties (Access = public)
        principalDirectionComputer
    end
    
    methods (Static) %(Access = {?Physical_Problem, ?Element_Elastic_Micro, ?obj})
        
        
        function obj = create(mesh,geometry,material,dof,problemData,interp)
            switch problemData.scale
                case 'MICRO'
                    switch problemData.pdim
                        case '2D'
                            obj = Element_Elastic_2D_Micro(mesh,geometry,material,dof,problemData,interp);
                        case '3D'
                            obj = Element_Elastic_3D_Micro(mesh,geometry,material,dof,problemData,interp);
                    end
                case 'MACRO'
                    switch problemData.pdim
                        case '2D'
                            obj = Element_Elastic_2D(mesh,geometry,material,dof,problemData,interp);
                        case '3D'
                            obj = Element_Elastic_3D(mesh,geometry,material,dof,problemData,interp);
                    end
            end
            obj.createPrincipalDirection(problemData.pdim);
        end
    end
    
    
    methods %(Access = {?Physical_Problem, ?Element_Elastic_Micro, ?Element})
        function compute(obj,mesh,geometry,material,dof,problemData,interp)
            obj.interpolation_u = interp{1};
            obj.initElement(geometry,mesh,material,dof,problemData.scale,interp)
            obj.nfields=1;
            
            obj.initialize_dvolum()
            Bmat   = obj.computeBmat();
            dvolum = obj.geometry.dvolu;
            
            ngaus   = obj.quadrature.ngaus;
            dimen   = obj.computeDim(ngaus);
            connect = mesh.connec;%obj.interp{1}.T;
            
            
            obj.dim = dimen;
            obj.connec = connect;
            %            obj.DeltaC = ContitutiveTensorIncrement();
            
            obj.K_generator = StiffnessMatrixGenerator(connect,Bmat,dvolum,dimen);
            
            obj.Bmatrix = obj.computeB_InMatrixForm();
            
            
            obj.StiffnessMatrix = KGeneratorWithfullStoredB(obj.dim,obj.connec,obj.Bmatrix,dvolum);
            
            
        end
        
        function dim = computeDim(obj,ngaus)
            dim                = DimensionVariables();
            dim.nnode          = obj.nnode;
            dim.nunkn          = obj.dof.nunkn;
            dim.nstre          = obj.nstre;
            dim.ndof           = obj.dof.ndof;
            dim.nelem          = obj.nelem;
            dim.ndofPerElement = dim.nnode*dim.nunkn;
            dim.ngaus          = ngaus;
            dim.nentries       = dim.nelem*(dim.ndofPerElement)^2;
        end
        
        function initialize_dvolum(obj)
            obj.computeQuadrature();
            obj.computeGeometry();
        end
        
        function computeQuadrature(obj)
            obj.quadrature.computeQuadrature('LINEAR');
        end
        
        function computeGeometry(obj)
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
        end
        
        function Kred = computeLHS(obj)
            obj.K = obj.computeStiffnessMatrix();
            Kred = obj.bcApplier.fullToReducedMatrix(obj.K);
        end
        
        function fext_red = computeRHS(obj)
            Fext = obj.computeExternalForces();
            R = obj.compute_imposed_displacement_force(obj.K);
            obj.fext = Fext + R;
            fext_red = obj.bcApplier.fullToReducedVector(obj.fext);
            obj.fextRed = fext_red;
        end
        
        function [K] = computeStiffnessMatrix(obj)
            [K] = obj.computeStiffnessMatrixSYM();
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
                                v = squeeze(Bmat(istre,iv,:).*Cmat(istre,jstre,:,igaus).*Bmat(jstre,jv,:));
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
            %  obj.K_generator.generate(obj.material.C);
            %  K = obj.K_generator.getStiffMatrix();
            
            
            obj.StiffnessMatrix.compute(obj.material.C);
            K = obj.StiffnessMatrix.K;
        end
        
        
        function Bmatrix = computeB_InMatrixForm(obj)
            
            ndofPerElement = obj.nnode*obj.dof.nunkn;
            Bfull = zeros(obj.quadrature.ngaus,obj.nstre,ndofPerElement,obj.nelem);
            Bmatrix = zeros(obj.quadrature.ngaus*obj.nelem*obj.nstre,ndofPerElement);
            for igaus = 1:obj.quadrature.ngaus
                unitaryIndex = false(obj.quadrature.ngaus*obj.nstre,1);
                pos = obj.nstre*(igaus-1) + 1 : obj.nstre*(igaus) ;
                unitaryIndex(pos) = true;
                
                Index = repmat(unitaryIndex,obj.nelem,1);
                Bfull(igaus,:,:,:) = obj.computeB(igaus);
                Bshif = reshape(permute(obj.computeB(igaus),[1 3 2]),obj.nelem*obj.nstre,ndofPerElement);
                Bmatrix(Index,:) = Bshif;
            end
        end
        
        function variables = computeVars(obj,uL)
            variables.d_u = obj.compute_displacements(uL);
            variables.fext = obj.fext;
            variables.strain = obj.computeStrain(variables.d_u,obj.dof.in_elem{1});
            variables.stress = obj.computeStress(variables.strain,obj.material.C,obj.quadrature.ngaus,obj.nstre);
            variables = obj.permuteStressStrain(variables);
            if ~isequal(class(obj),'Element_Elastic_3D')
                [dir,s] = obj.computePrincipalStressDirection(variables.stress);
                variables.principalDirections = dir;
                variables.principalStress     = s;
            end
        end
        
        
    end
    
    methods (Access = private)
        
        function createPrincipalDirection(obj, pdim)
            s.eigenValueComputer.type = 'PRECOMPUTED';
            s.type = pdim;
            p = PrincipalDirectionComputer.create(s);
            obj.principalDirectionComputer = p;
        end
        
        function [dir,str] = computePrincipalStressDirection(obj,tensor)
            obj.principalDirectionComputer.compute(tensor);
            dir = obj.principalDirectionComputer.direction;
            str = obj.principalDirectionComputer.principalStress;
        end
        
    end
    
    
    
    methods(Access = protected)
        function FextSuperficial = computeSuperficialFext(obj)
            FextSuperficial = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj)
            FextVolumetric = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
        
        function u = compute_displacements(obj,usol)
            u = obj.bcApplier.reducedToFullVector(usol);
        end
        
    end
    
    methods(Static)
        function variables = permuteStressStrain(variables)
            variables.strain = permute(variables.strain, [3 1 2]);
            variables.stress = permute(variables.stress, [3 1 2]);
        end
    end
    
    methods (Abstract, Access = protected)
        computeStrain(obj)
        computeBmat(obj)
    end
    
    methods (Abstract, Access = public)
        computeB(obj)
    end
    
    methods(Static, Access = protected)
        % Compute strains (e = B�u)
        
        % Compute stresses
        function stres = computeStress(strain,C,ngaus,nstre)
            stres = zeros(size(strain));
            for igaus = 1:ngaus
                for istre=1:nstre
                    for jstre=1:nstre
                        stres(istre,:,igaus) = stres(istre,:,igaus) + squeeze(C(istre,jstre,:,igaus))'.*strain(jstre,:,igaus);
                    end
                end
            end
        end
    end
end
