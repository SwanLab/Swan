classdef Element_Hyperelastic < Element
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
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
        end
        
        function [r,dr,K] = computeResidual(obj,uL)
            % *************************************************************
            % Compute
            % - residual: r = K(u)u - F
            % - residual derivative: dr = Ktan
            % *************************************************************
            
            % Compute tangent matrix
            K = obj.computeTangentMatrix();
            
            % Assemble
            K = obj.AssembleMatrix(K);
            
            % Assemble u and Fext
            u = zeros(obj.dof.ndof,1);
            u(obj.dof.vL) = uL;
            
            if ~isempty(obj.dof.vR)
                u(obj.dof.vR) = obj.bc.fixnodes(:,3);
                fext = obj.Fext(obj.dof.vL)-K(obj.dof.vL,obj.dof.vR)*u(obj.dof.vR); % fext + reac
            else
                fext = obj.Fext(obj.dof.vL);
            end
            
            %             % Compute internal forces
            %             fint = obj.computeInternal(sigma);
            %             fint = obj.AssembleVector(fint);
            %             fint = fint(obj.dof.vL);
            %
            %             r  = fint - fext;
            %             % ********************
            
            r = 0;
            dr = K(obj.dof.vL, obj.dof.vL);
            
            
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
%             fint = fint(obj.dof.vL);
        end
        
        
        function [K,sigma] = computeTangentMatrix(obj)
            % ctens   = obj.material.C;
            
            [ctens,sigma] = obj.material.computeCtens(obj.coord);
%             ctens = m.C;
%             sigma = m.sigma;
            
            

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
        
        function obj = updateCoord(obj,u)
            % Update coordinates
            coord0 = obj.coord;
            coord  = reshape(coord0(:,1:2)',[],1);
            coord(obj.dof.vL) = coord(obj.dof.vL) + u;
            coord0(:,1:2) = reshape(coord,2,[])';
            obj.coord = coord0;
            
            
        end
        
        function obj = updateCartd(obj,pdim)
            [obj.cartd,obj.dvolu] = obj.geometry.computeCartd(obj.coord,obj.nelem,pdim);
        end
        
    end
    
    
    methods(Access = protected)
        function FextSuperficial = computeSuperficialFext(obj,bc)
            FextSuperficial = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
    end
    
end

