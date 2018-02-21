classdef Element_Hyperelastic < Element
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fext
        cartd0
        cartd
        dvolu
    end
    
    methods (Access = {?Physical_Problem, ?Element})
        % Constructor
        function obj = Element_Hyperelastic(geometry)
            obj.cartd0 = geometry.cartd;
            
            % These variables will be updated
            obj.cartd = geometry.cartd;
            obj.dvolu = geometry.dvolu;
            
            % Nonlinear parameters
            obj.nincr = 85;
            obj.cload = 0;
        end
        
        function [r,dr,K] = computeResidual(obj,uL)
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
            
%             % Assemble u and Fext
%             u = zeros(obj.dof.ndof,1);
%             u(obj.dof.free) = uL;
%             
%             if ~isempty(obj.dof.vR)
%                 u(obj.dof.vR) = obj.bc.fixnodes(:,3);
%                 fext = obj.cload(obj.dof.vL)-K(obj.dof.vL,obj.dof.vR)*u(obj.dof.vR); % fext + reac
%             else
%                 fext = obj.cload(obj.dof.free);
%             end

            R = obj.compute_imposed_displacemet_force(K);
            obj.fext = obj.cload + R;

            fint_red = obj.computeInternal();
            
            r = fint_red - obj.fext(obj.dof.free);
            dr = K(obj.dof.free, obj.dof.free);
        end
        

        function fint = computeInternal(obj)
        % Fdef          --> cartd0 = obj.geometry.cartd
        % fint,ksigma   --> cartd  = obj.cartd
            sigma = obj.material.updateSigma(obj.coord,obj.cartd0);
            for a = 1:obj.nnode
                for i = 1:obj.geometry.ndime
                    iL = obj.geometry.ndime*(a-1) + i;
                    t = zeros(length(obj.dof.ndof),1,obj.nelem);
                    for j = 1:obj.geometry.ndime
                        t(j,:) = sigma(i,j,:).*obj.cartd(j,a,:);
                    end
                    t = squeeze(sum(t)).*squeeze(obj.dvolu);
                    t = permute(t,[3 2 1]);
                    T(iL,1,:) = t;
                end
            end
            fint = T;
            fint = obj.AssembleVector(fint);
            fint = fint(obj.dof.free);
        end
        
        
        function [K,sigma] = computeTangentMatrix(obj)
            % Compute ctens & sigma
            [ctens,sigma] = obj.material.computeCtens(obj.coord);

            % Compute tangent components
            kconst  = obj.computeConstitutive(ctens);
            ksigma  = obj.computeGeometric(sigma);
            K = kconst + ksigma;
        end
        
        
        function kconst = computeConstitutive(obj,ctens)
            % Initialization
            kconst = zeros(obj.nnode*obj.nunkn,obj.nnode*obj.nunkn,obj.nelem);
            
            % Constitutive component
            for a = 1:obj.nnode
                for b = 1:obj.nnode
                    for i = 1:obj.geometry.ndime
                        for j = 1:obj.geometry.ndime
                            iL = obj.geometry.ndime*(a-1) + i;
                            jL = obj.geometry.ndime*(b-1) + j;
                            for k = 1:obj.geometry.ndime
                                for l = 1:obj.geometry.ndime
                                    kconst(iL,jL,:) = kconst(iL,jL,:) + permute(squeeze(obj.cartd(k,a,:)).*squeeze(ctens(i,k,j,l,:)).*squeeze(obj.cartd(l,b,:)).*squeeze(obj.dvolu),[2 3 1]);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function ksigma = computeGeometric(obj,sigma)
            % Initialization
            ksigma = zeros(obj.nnode*obj.nunkn,obj.nnode*obj.nunkn,obj.nelem);
            
            % Vectorization identity matrix
            dk = eye(3);
            dk = repmat(dk,[1 1 obj.nelem]);
            
            % Initial stress component (geometric)
            for a = 1:obj.nnode
                for b = 1:obj.nnode
                    for i = 1:obj.geometry.ndime
                        for j = 1:obj.geometry.ndime
                            iL = obj.geometry.ndime*(a-1) + i;
                            jL = obj.geometry.ndime*(b-1) + j;
                            for k = 1:obj.geometry.ndime
                                for l = 1:obj.geometry.ndime
                                    ksigma(iL,jL,:) = ksigma(iL,jL,:) + permute(squeeze(obj.cartd(k,a,:).*sigma(k,l,:).*obj.cartd(l,b,:).*dk(i,j,:)).*squeeze(obj.dvolu),[2 3 1]);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % 2D
        function obj = updateCoord(obj,u)
            % Update coordinates
            coord0 = obj.coord;
            coord  = reshape(coord0(:,1:2)',[],1);
            coord(obj.dof.free) = coord(obj.dof.free) + u;
            coord0(:,1:2) = reshape(coord,2,[])';
            obj.coord = coord0;            
        end
        
        function obj = updateCartd(obj)
            [obj.cartd,obj.dvolu] = obj.geometry.computeCartd(obj.coord,obj.nelem,obj.pdim);
        end
        
        function variables = computeVars(obj,uL)
            variables = obj.computeDispStressStrain(uL);
        end
        
    end
    
    
    methods(Access = protected)  % Only the child sees the function
        function variables = computeDispStressStrain(obj,uL)
            variables.d_u = obj.compute_displacements(uL);
            variables.fext= obj.fext;
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

