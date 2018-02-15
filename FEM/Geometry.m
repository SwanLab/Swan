classdef Geometry
    %Geometry Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        cartd
        shape
        dvolu
        weigp
        ndime
        nnode
        ngaus
        djacb
        type
        deriv
        connec
    end
    properties (GetAccess = ?Postprocess, SetAccess = private)
        posgp
    end
    
    methods (Access = {?Physical_Problem, ?Element_Hyperelastic})
        function obj = Geometry(mesh)
            obj.nnode = size(mesh.connec,2);
            switch mesh.geometryType
                case 'TRIANGLE'
                    switch obj.nnode
                        case 3
                            geometryObject = Triangle_Linear;
                        case 6
                            geometryObject = Triangle_Quadratic;
                        otherwise
                            error('Invalid nnode for element TRIANGLE.');
                    end
                case 'Triangle_Linear_Mass'

                    geometryObject = Triangle_Linear_Mass;
                case 'Quad_Mass'
                    geometryObject = Quad_Mass;

                case 'QUAD'
                    switch obj.nnode
                        case 4
                            geometryObject = Quadrilateral_Bilinear;
                        case 8
                            geometryObject = Quadrilateral_Serendipity;
                        otherwise
                            error('Invalid nnode for element QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                    geometryObject = Tetrahedra;
                case 'HEXAHEDRA'
                    geometryObject = Hexahedra;
                otherwise
                    error('Invalid mesh type.')
            end
            obj.type  = geometryObject.type;
            obj.posgp = geometryObject.posgp;
            obj.weigp = geometryObject.weigp;
            obj.ndime = geometryObject.ndime;
            obj.ngaus = geometryObject.ngaus;
            obj.shape = geometryObject.shape;
            obj.deriv = geometryObject.deriv;
            obj.connec= mesh.connec;
            
            [cartd,dvolu,djacb] = obj.computeCartd(mesh.coord,mesh.nelem,mesh.pdim);
            
            obj.cartd = cartd;
            obj.dvolu = dvolu;
            obj.djacb = djacb;
        end
        
        function [cartd,dvolu,djacb] = computeCartd(obj,coord,nelem,pdim)
            for i = 1:obj.ndime
                a = coord(:,i);
                elcoord(:,i,:) = a(permute(obj.connec',[1,3,2]));
            end
            % Gauss loop
            for igauss = 1:obj.ngaus
                for i = 1:obj.ndime
                    % Jacobian
                    deriv_perm = permute(obj.deriv(i,:,igauss),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,obj.ndime,nelem);
                    jacobian(i,:,:) = sum(deriv_perm_large.*elcoord,1);
                end                
                % !! SWITCH EXECPCIÓ? !!
                switch pdim
                    case '2D'
                        [invJ,detJ] = multinverse2x2(jacobian);
                    case '3D'
                        [invJ,detJ] = multinverse3x3(jacobian);
                end
                
                for i = 1:obj.ndime
                    % Cartesian Derivatives
                    deriv_perm = permute(invJ(i,:,:),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,obj.nnode,1) .*repmat(obj.deriv(:,:,igauss),1,1,nelem);
                    cartd(i,:,:,igauss) = sum(deriv_perm_large);
                end
                dvolu(:,igauss) = obj.weigp(igauss)*detJ;
                djacb(:,igauss) = detJ;
            end 
        end
        
    end
    
end
% cartd
% dvolu
% djacb

