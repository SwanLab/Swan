classdef Geometry<handle
    %Geometry Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = public)
        type
        cartd
        cart_pos_gp
        dvolu
        djacob
        interpolation
        quadrature_order
        nfields=1;
    end
    methods %(Access = {?Physical_Problem,?Element_DiffReact}) % !! Element_DiffReact -> Chapusilla !!
        function obj = Geometry(mesh,order)
            obj.type = mesh.geometryType;
            obj.interpolation = Interpolation.create(mesh,order);
        end
        function computeGeometry(obj,quadrature,interp_variable)                
            if ~strcmp(obj.quadrature_order,quadrature.order)
                obj.compute(quadrature,interp_variable)
                obj.quadrature_order = quadrature.order;
            end
        end        
        function compute(obj,quadrature,interp_variable)
            ndime=interp_variable.ndime;
            nnode=interp_variable.nnode;
            ngaus=quadrature.ngaus;
            nelem=interp_variable.nelem;
            gp_position = zeros(ndime,ngaus,nelem);
            obj.dvolu = zeros(nelem,ngaus);
            obj.djacob = zeros(nelem,ngaus);
            obj.cartd=zeros(ndime,nnode,nelem,ngaus);
            for igauss=1:ngaus                
                for inode = 1:nnode
                    for idime=1:ndime
                        x = interp_variable.shape(inode,igauss).*squeeze(interp_variable.xpoints(interp_variable.T(:,inode),idime));
                        gp_position(idime,igauss,:) = squeeze(gp_position(idime,igauss,:))+x;
                    end
                end
            end

            obj.cart_pos_gp = gp_position;
            for i = 1:ndime
                a = interp_variable.xpoints(:,i);
                elcoord(:,i,:) = a(permute(interp_variable.T',[1,3,2]));
            end
            
            % Gauss loop
            for igauss = 1:ngaus
                for i = 1:ndime
                    % Jacobian
                    deriv_perm = permute(interp_variable.deriv(i,:,igauss),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,ndime,nelem);
                    jacobian(i,:,:) = sum(deriv_perm_large.*elcoord,1);
                end
                
                [invJ,detJ]=obj.inverseElementalMatrix(ndime,jacobian);
                
                for i = 1:ndime
                    % Cartesian Derivatives
                    deriv_perm = permute(invJ(i,:,:),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,interp_variable.nnode,1) .*repmat(interp_variable.deriv(:,:,igauss),1,1,nelem);
                    obj.cartd(i,:,:,igauss) = sum(deriv_perm_large);
                end
                obj.dvolu(:,igauss) = quadrature.weigp(igauss)*detJ;
                obj.djacob(:,igauss)= detJ;
            end
        end

    end
    
    methods (Static)
        function [inverse,determinant]=inverseElementalMatrix(ndime,A)
            switch ndime
                case 1
                    inverse = 1/A;
                    determinant = A;
                case 2
                    [inverse,determinant] = multinverse2x2(A);
                case 3
                    [inverse,determinant] = multinverse3x3(A);
            end
        end
    end
end