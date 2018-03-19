classdef Element_NeoHookean < Element_Hyperelastic
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % Constructor
        function obj = Element_NeoHookean(geometry)
            obj@Element_Hyperelastic(geometry); % when the child is created before the father
        end
        
        %% Compute Lagrangian Elasticity Tensor
        function Ctens = compute_Lagrangian_Elasticity(obj,Crcg,Fjacb)
            Ctens = zeros(3,3,3,3,obj.nelem);
            
            % Inverse of Crcg
            Cinv = multinverse3x3(Crcg);
            
            % Lame parameters
            mu      = obj.material.mu;
            lambda  = obj.material.lambda;
            
            % Material or Lagrangian Elasticity Tensor
            for I = 1:3
                for J = 1:3
                    for K = 1:3
                        for L = 1:3
                            Ctens(I,J,K,L,:) = squeeze(lambda*Cinv(I,J,:).*Cinv(K,L,:)) + 2*(mu - lambda*log(Fjacb)).*squeeze(1/2*(Cinv(I,K,:).*Cinv(J,L,:) + Cinv(I,L,:).*Cinv(J,K,:)));
                        end
                    end
                end
            end
        end
        
        %% Compute Second Piola-Kirchhoff Stress Tensor
        function S = compute_Second_PK_Stress(obj,Crcg,Fjacb)
            S = zeros(3,3,obj.nelem);
            
            % Inverse of Crcg
            Cinv = multinverse3x3(Crcg);
            
            % Identity matrix
            Imat = repmat(eye(3),[1 1 obj.nelem]);
            
            % Lame parameters
            mu      = obj.material.mu;
            lambda  = obj.material.lambda;
            
            for i = 1:3
                for j = 1:3
                    S(i,j,:) = mu*squeeze(Imat(i,j,:)-Cinv(i,j,:)) + lambda*squeeze(log(Fjacb)).*squeeze(Cinv(i,j,:));
                end
            end
        end
    end
    
end

