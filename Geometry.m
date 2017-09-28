classdef Geometry
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cartDeriv
        area
        ndime
        nnode % preguntar
        ngaus % preguntar
    end
    
    methods
        function obj = Geometry(mesh)
            switch mesh.geometryType
                case 'Triangle'
                    geometryObject = Triangle_Linear();
                case 'Tetrahedra'
                    geometryObject = Tetrahedra();
            end
            obj.ndime = geometryObject.ndime;
            obj.nnode = geometryObject.nnode;
            obj.ngaus=geometryObject.ngaus;
            for i = 1:obj.ndime
                a = mesh.coord(:,i);
                elcoord(:,i,:) = a(permute(mesh.connec',[1,3,2]));
            end
            for i = 1:obj.ndime
                % Jacobian
                deriv_perm = permute(geometryObject.deriv(i,:,:),[2,1,3]);
                deriv_perm_large = repmat(deriv_perm,1,obj.ndime,mesh.nelem);
                jacobian(i,:,:) = sum(deriv_perm_large.*elcoord,1);
            end
            switch obj.ndime
                case 2
                    [invJ,detJ] = multinverse2x2(jacobian);
                case 3
                    [invJ,detJ] = multinverse3x3(jacobian);
            end
            for i = 1:obj.ndime
                % Cartesian Derivatives
                deriv_perm = permute(invJ(i,:,:),[2,1,3]);
                deriv_perm_large = repmat(deriv_perm,1,geometryObject.nnode,1) .*repmat(geometryObject.deriv,1,1,mesh.nelem);
                obj.cartDeriv(i,:,:) = sum(deriv_perm_large);
            end
            % Area
            obj.area = geometryObject.weigp*detJ;
        end
    end
    
end

