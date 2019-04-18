classdef Isotropic2dHyperElasticMaterial < ElasticMaterial
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Element_Hyperelastic,?PhysicalVars_Elastic}, SetAccess = protected)
        mup
        lambdap
        connec
        cartd
        nnode
        coord
        sigma
    end
    
    methods
        function obj = Isotropic2dHyperElasticMaterial(cParams)
            obj.nelem = cParams.nelem;
            obj.connec= cParams.connec;
            obj.cartd = cParams.cartd;
            obj.nnode = cParams.nnode;
%             obj.coord = coord(:,1:2);
            obj = computeCtens(obj,cParams.coord(:,1:2)); % obj.computeCtens
%             obj = obj.updateConstitutiveTensor(coord);
        end
        
        %% Compute Eulerian elasticity tensor
        function obj = computeCtens(obj,coord)
            
            % Coordinate's vectorization
            x = coord(obj.connec(:,:)',:)';
            x = reshape(x,[2,obj.nnode,obj.nelem]);
            
            % Jacobian
            [Fjacb,blcg] = obj.computeDefGradient(x,obj.cartd,obj.nnode,obj.nelem);
            
            % Effective Lame moduli
            obj.mup     = (obj.mu-obj.lambda*log(Fjacb))./Fjacb;
            obj.lambdap = obj.lambda./Fjacb;
            
            % Rearrange the dimensions
            obj.mu      = repmat(obj.mu,[1 1 2]);
            obj.mup     = permute(obj.mup,[2 3 1]);
            obj.lambda  = repmat(obj.lambda,[1 1 2]);
            obj.lambdap = permute(obj.lambdap,[2 3 1]);
            
            % Define constitutive tensor
            obj.C  = zeros(3,3,3,3,obj.nelem);
            
            % Define delta kronecker
            dk = eye(3);
            dk = repmat(dk,[1 1 obj.nelem]);
            
            % Cauchy stress tensor for a compressible neo-Hookean material 
            obj = obj.computeNeoHookean(blcg,Fjacb);
            
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
        function obj = computeNeoHookean(obj,blcg,Fjacb)
            I = repmat(eye(3),[1 1 obj.nelem]);
            Fjacb = permute(Fjacb,[2 3 1]);
            Fjacb = repmat(Fjacb,3,3);
            obj.mu = repmat(obj.mu,3,3);
            obj.lambda = repmat(obj.lambda,3,3);
            obj.sigma = obj.mu./Fjacb.*(blcg-I) + obj.lambda./Fjacb.*(log(Fjacb)).*I;
        end
        
        
        %% Update material properties
        function obj = updateConstitutiveTensor(obj,coord)
            obj = obj.computeCtens(coord);
        end
    end
    
    methods (Static)
        %% Compute deformation gradient tensor
        function [Fjacb,blcg] = computeDefGradient(x,cartd0,nnode,nelem)
            F = repmat(eye(3),[1 1 nelem]);
            for i = 1:2 % 2D
                f = zeros(2,nnode,nelem);
                for j = 1:2
                    for a = 1:nnode
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