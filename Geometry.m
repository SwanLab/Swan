classdef Geometry
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cartDeriv
        area
        ndime
        nnode
    end
    
    
    methods
        function obj = Geometry(mesh)
            triangleLinear = Triangle_Linear();
            obj.ndime = triangleLinear.ndime;
            obj.nnode = triangleLinear.nnode;
            for i = 1:obj.ndime
                a = mesh.coord(:,i);
                elcoord(:,i,:) = a(permute(mesh.connec',[1,3,2]));
            end
            for i = 1:obj.ndime
                % Jacobian
                deriv_perm = permute(triangleLinear.deriv(i,:,:),[2,1,3]);
                deriv_perm_large = repmat(deriv_perm,1,obj.ndime,mesh.nelem);
                jacobian(i,:,:) = sum(deriv_perm_large.*elcoord,1);
            end
            
            [invJ,detJ] = multinverse2x2(jacobian);
            
            for i = 1:obj.ndime
                % Cartesian Derivatives
                deriv_perm = permute(invJ(i,:,:),[2,1,3]);
                deriv_perm_large = repmat(deriv_perm,1,triangleLinear.nnode,1) .*repmat(triangleLinear.deriv,1,1,mesh.nelem);
                obj.cartDeriv(i,:,:) = sum(deriv_perm_large);
            end
            % Area
            obj.area = triangleLinear.weigp*detJ;
        end
    end
    
end

