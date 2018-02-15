classdef Material_Hyperelastic_2D < Material_Elastic
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Element_Hyperelastic,?PhysicalVars_Elastic, ?Physical_Problem}, SetAccess = ?Physical_Problem)
        mup
        lambdap
        connec
        cartd
        nnode
        coord
        sigma
    end
    
    methods
        function obj = Material_Hyperelastic_2D(nelem,connec,cartd,nnode,coord)
            obj@Material_Elastic(nelem);
            obj.connec= connec;
            obj.cartd = cartd;
            obj.nnode = nnode;
            %             obj = computeCtens(obj,coord(:,1:2)); % obj.computeCtens
            
            obj.coord = coord;
            
            % Jacobian
            [Fjacb,~] = obj.computeDefGradient(coord,obj.cartd);
            
            % Effective Lame moduli
            obj.mup     = (obj.mu-obj.lambda*log(Fjacb))./Fjacb;
            obj.lambdap = obj.lambda./Fjacb;
            
            % Vectorize
            obj.mu      = repmat(obj.mu,[1 1 nelem]);
            obj.lambda  = repmat(obj.lambda,[1 1 nelem]);
            
            % Rearrange the dimensions
            obj.mup     = permute(obj.mup,[2 3 1]);
            obj.lambdap = permute(obj.lambdap,[2 3 1]);
        end
        
        %% Compute Eulerian elasticity tensor
        function obj = computeCtens(obj,coord)
            
            % Jacobian
            [Fjacb,blcg] = obj.computeDefGradient(coord,obj.cartd);
            
            % Cauchy stress tensor for a compressible neo-Hookean material
            obj.sigma = obj.computeNeoHookean(blcg,Fjacb);
            
            % Define constitutive tensor
            obj.C  = zeros(3,3,3,3,obj.nelem);
            
            % Define delta kronecker
            dk = eye(3);
            dk = repmat(dk,[1 1 obj.nelem]);
            
            % Compute c
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l = 1:3
                            inc = obj.lambdap.*dk(i,j,:).*dk(k,l,:) + obj.mup.*(dk(i,k,:).*dk(j,l,:) + dk(i,l,:).*dk(j,k,:));
                            obj.C(i,j,k,l,:) = squeeze(obj.C(i,j,k,l,:)) + squeeze(inc);
                        end
                    end
                end
            end
        end
        
        %% Evaluate Cauchy stress tensor for a compressible neo-Hookean material
        function sigma = computeNeoHookean(obj,blcg,Fjacb)
            I = repmat(eye(3),[1 1 obj.nelem]);
            Fjacb = permute(Fjacb,[2 3 1]);
            Fjacb = repmat(Fjacb,3,3);
            obj.mu = repmat(obj.mu,3,3);
            obj.lambda = repmat(obj.lambda,3,3);
            sigma = obj.mu./Fjacb.*(blcg-I) + obj.lambda./Fjacb.*(log(Fjacb)).*I;
        end
        
        
        %% Update material properties
        function sigma = updateSigma(obj,coord,cartd0)
            [Fjacb,blcg] = obj.computeDefGradient(coord,cartd0);
            sigma = obj.computeNeoHookean(blcg,Fjacb);
        end
        
        %% Compute deformation gradient tensor
        % Compute F, b & sigma
        function [Fjacb,blcg] = computeDefGradient(obj,coord,cartd0)
            % 2D
            coord = coord(:,1:2);
            
            % Coordinate's vectorization
            x = coord(obj.connec(:,:)',:)';
            x = reshape(x,[2,obj.nnode,obj.nelem]);
            
            % Deformation gradient tensor
            F = repmat(eye(3),[1 1 obj.nelem]);
            
            for i = 1:2 % 2D
                f = zeros(2,obj.nnode,obj.nelem);
                for j = 1:2
                    for a = 1:obj.nnode
                        inc = x(i,a,:).*cartd0(j,a,:);
                        f(j,a,:) = f(j,a,:) + inc;
                    end
                    F(i,j,:) = sum(f(j,:,:));
                end
            end
            
            % Determinant vectorization
            [~,Fjacb] = multinverse3x3(F);
            
            % Left-Cauchy deformation tensor
            blcg = F.*permute(F,[2 1 3]);
        end
    end
    
end