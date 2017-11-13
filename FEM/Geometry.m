classdef Geometry
    %Geometry Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        cartDeriv
        area
        ndime
        nnode % preguntar
        ngaus % preguntar
    end
    
    methods (Access = ?Physical_Problem)
        function obj = Geometry(mesh)
            switch mesh.geometryType
                case 'TRIANGLE'
                    geometryObject = Triangle_Linear();
                case 'QUAD'
                    geometryObject = Quadrilateral();
                case 'TETRAHEDRA'
                    geometryObject = Tetrahedra();
                case 'HEXAHEDRA'
                    geometryObject = Hexahedra();
            end
            
            obj.ndime = geometryObject.ndime;
            obj.nnode = geometryObject.nnode;
            obj.ngaus = geometryObject.ngaus;
            for i = 1:obj.ndime
                a = mesh.coord(:,i);
                elcoord(:,i,:) = a(permute(mesh.connec',[1,3,2]));
            end
            %gauss loop
            for igauss=1:obj.ngaus
                for i = 1:obj.ndime
                    % Jacobian
                    deriv_perm = permute(geometryObject.deriv(i,:,igauss),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,obj.ndime,mesh.nelem);
                    jacobian(i,:,:) = sum(deriv_perm_large.*elcoord,1);
                end
                switch obj.ndime
                    case 1
                        invJ = 1./jacobian;
                        detJ = jacobian;
                    case 2
                        [invJ,detJ] = multinverse2x2(jacobian);
                    case 3
                        [invJ,detJ] = multinverse3x3(jacobian);
                end
                for i = 1:obj.ndime
                    % Cartesian Derivatives
                    deriv_perm = permute(invJ(i,:,:),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,geometryObject.nnode,1) .*repmat(geometryObject.deriv(:,:,igauss),1,1,mesh.nelem);
                    obj.cartDeriv(i,:,:,igauss) = sum(deriv_perm_large);
                end
                obj.area(:,igauss) = geometryObject.weigp(igauss)*detJ;
            end
            % Area
            
            
        end
    end
    
end

