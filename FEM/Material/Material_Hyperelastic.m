classdef Material_Hyperelastic < Material_Elastic
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Element_Hyperelastic, ?Material_Hyperelastic_2D, ?Material_Hyperelastic_3D, ?PhysicalVars_Elastic, ?Physical_Problem}, SetAccess = ?Physical_Problem)
        connec
        cartd0
        nnode
        coord
    end
    
    methods
        function obj = Material_Hyperelastic(nelem,connec,cartd,nnode,coord)
            obj@Material_Elastic(nelem);
            obj.connec= connec;
            obj.cartd0 = cartd;
            obj.nnode = nnode;
            obj.coord = coord;
        end
        
        %% Compute Eulerian elasticity tensor
        function [ctens,sigma] = computeCtens(obj,coord,igaus)
            % Jacobian
            [F,Fjacb,Crcg] = obj.computeDefGradient(coord,obj.cartd0(:,:,:,igaus));
            
            % Neo-hookean -> 'neo'
            % St. Venant  -> 'stv'
            type = 'neo';
            
            % Material elasticity tensor & Second Piola-Kirchhoff stress
            [Ctens,S] = obj.computeElasticity(Fjacb,obj.lambda,obj.mu,Crcg,type);
            
            % Spatial elasticity tensor & Cauchy stress
            [ctens,sigma] = obj.pushForward(Ctens,S,F,Fjacb);
        end
        
        %% Compute Material elasticity & stress
        function [Ctens,S] = computeElasticity(obj,Fjacb,lambda,mu,Crcg,type)
            dk = eye(3);
            dk = repmat(dk,[1 1 obj.nelem]);
            
            Cinv = multinverse3x3(Crcg);
            Imat = repmat(eye(3),[1 1 obj.nelem]);
            
            Ctens = zeros(3,3,3,3,obj.nelem);
            S = zeros(3,3,obj.nelem);
            
            % Strain
            E = 1/2*(Crcg - Imat);
            
            switch type
                case 'neo'
                    
                    % Elasticity tensor - NEOHOOKEAN
                    for I = 1:3
                        for J = 1:3
                            for K = 1:3
                                for L = 1:3
                                    Ctens(I,J,K,L,:) = squeeze(lambda*Cinv(I,J,:).*Cinv(K,L,:)) + 2*(mu - lambda*log(Fjacb)).*squeeze(1/2*(Cinv(I,K,:).*Cinv(J,L,:) + Cinv(I,L,:).*Cinv(J,K,:)));
                                end
                            end
                        end
                    end
                    
                    for i = 1:3
                        for j = 1:3
                            S(i,j,:) = obj.mu*squeeze(Imat(i,j,:)-Cinv(i,j,:)) + obj.lambda*squeeze(log(Fjacb)).*squeeze(Cinv(i,j,:));
                        end
                    end
                    
                case 'stv'
                    
                    
                    % Elasticity tensor - ST.VENANT
                    for I = 1:3
                        for J = 1:3
                            for K = 1:3
                                for L = 1:3
                                    Ctens(I,J,K,L,:) = lambda*dk(I,J,:).*dk(K,L,:) + mu*(dk(I,K,:).*dk(J,L,:) + dk(I,L,:).*dk(J,K,:));
                                end
                            end
                        end
                    end
                    
                    for i = 1:3
                        for j = 1:3
                            for k = 1:3
                                S(i,j,:) = S(i,j,:) + obj.lambda*E(k,k,:)*dk(i,j);
                            end
                            S(i,j,:) = S(i,j,:) + 2*obj.mu*E(i,j,:);
                        end
                    end
            end
            
            
        end
        
        %% Push forward
        function [ctens,sigma] = pushForward(obj,Ctens,S,F,Fjacb)
            ctens = zeros(3,3,3,3,obj.nelem);
            sigma = zeros(3,3,obj.nelem);
            
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l = 1:3
                            sigma(i,l,:) = squeeze(sigma(i,l,:)) + 1./Fjacb.*squeeze(F(i,k,:)).*squeeze(S(k,j,:)).*squeeze(F(l,j,:));
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
                                            ctens(i,j,k,l,:) = squeeze(ctens(i,j,k,l,:)) + 1./Fjacb.*squeeze(F(i,I,:).*F(j,J,:).*F(k,K,:).*F(l,L,:)).*squeeze(Ctens(I,J,K,L,:));
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