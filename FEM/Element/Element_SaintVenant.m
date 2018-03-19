classdef Element_SaintVenant < Element_Hyperelastic
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % Constructor
        function obj = Element_SaintVenant(geometry)
            obj@Element_Hyperelastic(geometry); % when the child is created before the father
        end
        
        %% Compute Lagrangian Elasticity Tensor
        function Ctens = compute_Lagrangian_Elasticity(obj,~,~)
            Ctens = zeros(3,3,3,3,obj.nelem,obj.geometry.ngaus);
            % Kronecker delta
            dk = eye(3);
            dk = repmat(dk,[1 1 obj.nelem obj.geometry.ngaus]);
            
            % Lame parameters
            mu      = obj.material.mu;
            lambda  = obj.material.lambda;
            
            % Material or Lagrangian Elasticity Tensor
            for igaus = 1:obj.geometry.ngaus
                for I = 1:3
                    for J = 1:3
                        for K = 1:3
                            for L = 1:3
                                Ctens(I,J,K,L,:,igaus) = lambda*dk(I,J,:,igaus).*dk(K,L,:,igaus) + mu*(dk(I,K,:,igaus).*dk(J,L,:,igaus) + dk(I,L,:,igaus).*dk(J,K,:,igaus));
                            end
                        end
                    end
                end
            end
        end
        
        %% Compute Second Piola-Kirchhoff Stress Tensor
        function S = compute_Second_PK_Stress(obj,Crcg,~)
            S = zeros(3,3,obj.nelem,obj.geometry.ngaus);
            
            % Kronecker delta
            dk = eye(3);
            dk = repmat(dk,[1 1 obj.nelem]);
            
            % Lame parameters
            mu      = obj.material.mu;
            lambda  = obj.material.lambda;
            
            % Identity matrix
            Imat = repmat(eye(3),[1 1 obj.nelem]);
            
            for igaus = 1:obj.geometry.ngaus
                
                E = 1/2*(Crcg(:,:,:,igaus) - Imat); % Lagrangian Strain
                
                for i = 1:3
                    for j = 1:3
                        for k = 1:3
                            S(i,j,:,igaus) = S(i,j,:,igaus) + lambda*E(k,k,:)*dk(i,j);
                        end
                        S(i,j,:,igaus) = S(i,j,:,igaus) + 2*mu*E(i,j,:);
                    end
                end
            end
        end
    end
    
end

