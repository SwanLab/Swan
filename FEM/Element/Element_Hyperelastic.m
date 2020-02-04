classdef Element_Hyperelastic < Element
    %Element_Hyperelastic Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Access = public)
        function [r,dr] = computeResidual(obj,uL)
            % *************************************************************
            % Compute
            % - residual: r = K(u)u - F
            % - residual derivative: dr = Ktan
            % *************************************************************
            
            % Compute tangent matrix
            [K] = obj.computeTangentMatrix();
            
            % Assemble
            [K] = obj.AssembleMatrix(K);
            
            %Assemble u and Fext
            u = zeros(obj.dof.ndof,1);
            u(obj.dof.vL) = uL;
            
            if ~isempty(obj.dof.vR)
                u(obj.dof.vR) = obj.bc.fixnodes(:,3);
                fext = obj.Fext(obj.dof.vL)-K(obj.dof.vL,obj.dof.vR)*u(obj.dof.vR); % fext + reac
            else
                fext = obj.Fext(obj.dof.vL);
            end
            
            % Compute internal forces
            fint = obj.computeInternal();
            fint = obj.AssembleVector(fint);
            
            dr = K;            
            nincr = 100;
            fincr = fext/nincr;
            for incrm = 1:nincr
                r = fint - fext;
            end
        end
        
        
        function K = computeTangentMatrix(obj)
            ctens   = obj.material.C;
            kconst  = obj.computeConstitutive(ctens);
            ksigma  = obj.computeGeometric();
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
                                    kconst(iL,jL,:) = squeeze(kconst(iL,jL,:)) + squeeze(obj.material.cartd(k,a,:).*ctens(i,k,j,l).*obj.material.cartd(l,b,:)).*squeeze(obj.geometry.dvolu);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function ksigma = computeGeometric(obj)
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
                                    ksigma(iL,jL,:) = squeeze(ksigma(iL,jL,:)) + squeeze(obj.material.cartd(k,a,:).*obj.material.sigma(k,l,:).*obj.material.cartd(l,b).*dk(i,j,:)).*squeeze(obj.geometry.dvolu);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    methods(Access = protected)
        function FextSuperficial = computeSuperficialFext(obj,bc)
            FextSuperficial = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        function fint = computeInternal(obj)
            for a = 1:obj.nnode
                for i = 1:obj.geometry.ndime
                    iL = obj.geometry.ndime*(a-1) + i;
                    t = zeros(length(obj.dof.ndof),1,obj.nelem);
                    for j = 1:obj.geometry.ndime
                        t(j,:) = obj.material.sigma(i,j,:).*obj.material.cartd(j,a,:);
                    end
                    t = squeeze(sum(t)).*squeeze(obj.geometry.dvolu);
                    t = permute(t,[3 2 1]);
                    T(iL,1,:) = t;
                    
                end
            end
            fint = T;
        end
    end
    
end

