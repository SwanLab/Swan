classdef Element_Hyperelastic < Element
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fext
        cartd0
        cartd
        dvolu
        stress
    end
    
    methods (Access = {?Physical_Problem, ?Element})
        % Constructor
        function obj = Element_Hyperelastic(geometry)
            obj.cartd0 = geometry.cartd;
            
            % These variables will be updated
            obj.cartd = geometry.cartd;
            obj.dvolu = geometry.dvolu;
            
            % Nonlinear parameters
            obj.nincr = 80;
        end
        
        function [r,dr,K,fint] = computeResidual(obj,uL,lambda_i)
            % *************************************************************
            % Compute
            % - residual: r = K(u)u - F
            % - residual derivative: dr = Ktan
            % *************************************************************
            
            % Only for Element_Hyperelastic
            obj.updateCoord(uL);
            obj.updateCartd();
            
            % Compute tangent matrix
            K = obj.computeTangentMatrix();
            
            % Assemble
            K = obj.AssembleMatrix(K);
            
            R = obj.compute_imposed_displacemet_force(K);
            obj.fext = obj.Fext + R; % fext + reac

            % Compute internal forces
            [fint_red, fint] = obj.computeInternal();

            r = fint_red - lambda_i*obj.fext(obj.dof.free);
            dr = K(obj.dof.free, obj.dof.free);
        end
        

        function [fint_red, fint] = computeInternal(obj)
            T = zeros(length(obj.dof.in_elem(:,1)),1,obj.nelem);
            for igaus = 1:obj.geometry.ngaus
                for a = 1:obj.nnode
                    for i = 1:obj.geometry.ndime
                        iL = obj.geometry.ndime*(a-1) + i;
                        t = zeros(length(obj.dof.ndof),1,obj.nelem);
                        for j = 1:obj.geometry.ndime
                            t(j,:) = obj.stress(i,j,:,igaus).*obj.cartd(j,a,:,igaus);
                        end
                        t = squeeze(sum(t)).*squeeze(obj.dvolu(:,igaus));
                        t = permute(t,[3 2 1]);
                        Ti(iL,1,:) = t;
                    end
                end
                T = T + Ti;
            end
            fint = T;
            fint = obj.AssembleVector(fint);
            fint_red = fint(obj.dof.free);
        end
        
        
        function [K,sigma] = computeTangentMatrix(obj)
                % Compute ctens & sigma
                [ctens,sigma] = obj.compute_Hyperelasticity(obj.coord);
                obj.stress = sigma;

                % Compute tangent components
                kconst  = obj.computeConstitutive(ctens);
                ksigma  = obj.computeGeometric(sigma);
                K = kconst + ksigma;
        end
        
        % Constitutive component
        function kconst = computeConstitutive(obj,ctens)
            % Initialization
            kconst = zeros(obj.nnode*obj.nunkn,obj.nnode*obj.nunkn,obj.nelem);
            
            for igaus = 1:obj.geometry.ngaus
                for a = 1:obj.nnode
                    for b = 1:obj.nnode
                        for i = 1:obj.geometry.ndime
                            for j = 1:obj.geometry.ndime
                                iL = obj.geometry.ndime*(a-1) + i;
                                jL = obj.geometry.ndime*(b-1) + j;
                                for k = 1:obj.geometry.ndime
                                    for l = 1:obj.geometry.ndime
                                        %         kconst(iL,jL,:) = kconst(iL,jL,:) + permute(squeeze(obj.cartd(k,a,:)).*squeeze(ctens(i,k,j,l,:)).*squeeze(obj.cartd(l,b,:)).*squeeze(obj.dvolu),[2 3 1]);
                                        add1 = permute(kconst(iL,jL,:),[3 1 2]);
                                        add2 = permute(obj.cartd(k,a,:,igaus),[3 2 1]).*permute(ctens(i,k,j,l,:,igaus),[5 1 2 3 4]).*permute(obj.cartd(l,b,:,igaus),[3 2 1]).*permute(obj.dvolu(:,igaus),[1 3 2]);
                                        kconst(iL,jL,:) = bsxfun(@plus,add1,add2);
                                        %         kconst(iL,jL,:) = permute(kconst(iL,jL,:),[3 1 2]) + permute(obj.cartd(k,a,:),[3 2 1]).*permute(ctens(i,k,j,l,:),[5 1 2 3 4]).*permute(obj.cartd(l,b,:),[3 2 1]).*permute(obj.dvolu,[1 3 2]);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % Initial stress component (geometric)
        function ksigma = computeGeometric(obj,sigma)
            % Initialization
            ksigma = zeros(obj.nnode*obj.nunkn,obj.nnode*obj.nunkn,obj.nelem);
            
            % Vectorization identity matrix
            dk = eye(3); dk = repmat(dk,[1 1 obj.nelem]);
            
            for igaus = 1:obj.geometry.ngaus
                for a = 1:obj.nnode
                    for b = 1:obj.nnode
                        for i = 1:obj.geometry.ndime
                            for j = 1:obj.geometry.ndime
                                iL = obj.geometry.ndime*(a-1) + i;
                                jL = obj.geometry.ndime*(b-1) + j;
                                for k = 1:obj.geometry.ndime
                                    for l = 1:obj.geometry.ndime
                                        add1 = permute(ksigma(iL,jL,:),[3 1 2]);
                                        add2 = permute(obj.cartd(k,a,:,igaus).*sigma(k,l,:,igaus).*obj.cartd(l,b,:,igaus).*dk(i,j,:),[3 1 2]).*obj.dvolu(:,igaus);
                                        ksigma(iL,jL,:) = bsxfun(@plus,add1,add2);
                                        %                                     ksigma(iL,jL,:) = permute(ksigma(iL,jL,:),[3 1 2]) + permute(obj.cartd(k,a,:).*sigma(k,l,:).*obj.cartd(l,b,:).*dk(i,j,:),[3 1 2]).*obj.dvolu;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        
        %% Update Functions ***********************************************
        
        % Update Cartesian Derivatives
        function obj = updateCartd(obj)
            [obj.cartd,obj.dvolu] = obj.geometry.computeCartd(obj.coord,obj.nelem,obj.pdim);
        end
        
        % Update Position
        function obj = updateCoord(obj,u)
            coord0 = obj.coord;
            coord  = reshape(coord0(:,1:obj.geometry.ndime)',[],1);
            
            % Update
            coord(obj.dof.free) = coord(obj.dof.free) + u;
            
            % Correct format
            coord0(:,1:obj.geometry.ndime) = reshape(coord,obj.geometry.ndime,[])';
            obj.coord = coord0;            
        end
        
        %% Compute Variables **********************************************
        
        function variables = computeVars(obj,uL)
            variables = obj.computeDispStressStrain(uL);
        end
        
        %% Hyperelasticity ************************************************
        
        % Compute Eulerian Elasticity and main tensors
        function [ctens,sigma] = compute_Hyperelasticity(obj,coord)
            % Deformation Gradient tensor
            [F,Fjacb] = obj.compute_Deformation_Gradient(coord,obj.cartd0);
            
            % Right Cauchy Deformation tensor
            Crcg = obj.compute_Right_Cauchy_Deformation(F);
            
            % Lagrangian Elasticity tensor
            Ctens = obj.compute_Lagrangian_Elasticity(Crcg,Fjacb);
            
            % Second Piola-Kirchhoff Stress tensor
            S = obj.compute_Second_PK_Stress(Crcg,Fjacb);
            
            % Eulerian Elasticity & Cauchy Stress tensors
            [ctens,sigma] = obj.pushForward(Ctens,S,F,Fjacb);
        end
        
        % Compute deformation gradient tensor
        function [F,Fjacb] = compute_Deformation_Gradient(obj,coord,cartd0)
            coord = coord(:,1:obj.geometry.ndime);
            
            % Coordinate's vectorization
            x = coord(obj.geometry.connec(:,:)',:)';
            x = reshape(x,[obj.geometry.ndime,obj.nnode,obj.nelem]);
            
            % Deformation gradient tensor
            Fi = repmat(eye(3),[1 1 obj.nelem]);
            F  = zeros(3,3,obj.nelem);
            
            for igaus = 1:obj.geometry.ngaus
                for i = 1:obj.geometry.ndime
                    f = zeros(obj.geometry.ndime,obj.nnode,obj.nelem);
                    for j = 1:obj.geometry.ndime
                        for a = 1:obj.nnode
                            inc = x(i,a,:).*cartd0(j,a,:,igaus);
                            f(j,a,:) = f(j,a,:) + inc;
                        end
                        Fi(i,j,:) = sum(f(j,:,:));
                    end
                end
                F(:,:,:,igaus) = Fi;
                
             % Determinant vectorization (Jacobian)
            [~,Fjacbi] = multinverse3x3(Fi);
            Fjacb(:,igaus) = Fjacbi;
            end
        end
        
        % Compute Right-Cauchy Deformation Gradient
        function Crcg = compute_Right_Cauchy_Deformation(obj,F)
            Crcg = zeros(3,3,obj.nelem,obj.geometry.ngaus);
            for igaus = 1:obj.geometry.ngaus
                for i = 1:3
                    for j = 1:3
                        for k = 1:3
                            Crcg(i,j,:,igaus) = squeeze(Crcg(i,j,:,igaus)) + (squeeze(F(k,i,:,igaus))).*(squeeze(F(k,j,:,igaus)));
                        end
                    end
                end
            end
        end
        
        % Push forward
        function [ctens,sigma] = pushForward(obj,Ctens,S,F,Fjacb)
            ctens = zeros(3,3,3,3,obj.nelem,obj.geometry.ngaus);
            sigma = zeros(3,3,obj.nelem,obj.geometry.ngaus);
            
            for igaus = 1:obj.geometry.ngaus
                for i = 1:3
                    for j = 1:3
                        for k = 1:3
                            for l = 1:3
                                sigma(i,l,:,igaus) = permute(sigma(i,l,:,igaus),[3 2 1]) + 1./Fjacb(:,igaus).*permute(F(i,k,:,igaus).*S(k,j,:,igaus).*F(l,j,:,igaus),[3 2 1]);
                            end
                        end
                    end
                end
                
                for i = 1:3
                    for j = 1:3
                        for k = 1:3
                            for l = 1:3
                                for I = 1:3
                                    for J = 1:3
                                        for K = 1:3
                                            for L = 1:3
                                                add1 = permute(ctens(i,j,k,l,:,igaus),[5 1 2 3 4]);
                                                mul1 = 1./Fjacb(:,igaus).*permute(F(i,I,:,igaus).*F(j,J,:,igaus).*F(k,K,:,igaus).*F(l,L,:,igaus),[3 2 1]);
                                                mul2 = permute(Ctens(I,J,K,L,:,igaus),[5 1 2 3 4]);
                                                add2 = bsxfun(@times,mul1,mul2);
                                                ctens(i,j,k,l,:,igaus) = bsxfun(@plus,add1,add2);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    
    methods(Access = protected)  % Only the child sees the function
        function variables = computeDispStressStrain(obj,uL)
            variables.d_u = obj.compute_displacements(uL);
            % variables.strain
            variables.fext= obj.fext;
            variables.stress = obj.stress;
        end
        
        function u = compute_displacements(obj,usol)
            u = obj.reduced_vector_2_full_vector(usol,obj.dof);
        end
        
        function FextSuperficial = computeSuperficialFext(obj,bc)
            FextSuperficial = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
    end
    
    methods(Static)
        function Ared = full_matrix_2_reduced_matrix(A,dof)
            Ared = A(dof.free,dof.free);
        end
        
        function b_red = full_vector_2_reduced_vector(b,dof)
            b_red = b(dof.free);
        end
        
        function b = reduced_vector_2_full_vector(bfree,dof)
            b = zeros(dof.ndof,1);
            b(dof.free) = bfree;
            b(dof.dirichlet) = dof.dirichlet_values;
        end 
    end
    
end
    
