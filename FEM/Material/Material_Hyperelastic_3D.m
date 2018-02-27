classdef Material_Hyperelastic_3D < Material_Hyperelastic
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Element_Hyperelastic,?PhysicalVars_Elastic, ?Physical_Problem}, SetAccess = ?Physical_Problem)
    end
    
    methods
        function obj = Material_Hyperelastic_3D(nelem,connec,cartd,nnode,coord)
            obj@Material_Hyperelastic(nelem,connec,cartd,nnode,coord);
        end
        %% Compute deformation gradient tensor
        function [F,Fjacb,Crcg] = computeDefGradient(obj,coord,cartd0)
            
            % Coordinate's vectorization
            x = coord(obj.connec(:,:)',:)';
            x = reshape(x,[3,obj.nnode,obj.nelem]);
            
            % Deformation gradient tensor
            F = repmat(eye(3),[1 1 obj.nelem]);
            
            for i = 1:3 % 3D
                f = zeros(2,obj.nnode,obj.nelem);
                for j = 1:3
                    for a = 1:obj.nnode
                        inc = x(i,a,:).*cartd0(j,a,:);
                        f(j,a,:) = f(j,a,:) + inc;
                    end
                    F(i,j,:) = sum(f(j,:,:));
                end
            end
            
            % Determinant vectorization (Jacobian)
            [~,Fjacb] = multinverse3x3(F);
            
            % Right-Cauchy deformation tensor
            Crcg = zeros(3,3,obj.nelem);
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        Crcg(i,j,:) = squeeze(Crcg(i,j,:)) + (squeeze(F(k,i,:))).*(squeeze(F(k,j,:)));
                    end
                end
            end
        end
        
    end
end