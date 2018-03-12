classdef Geometry<handle
    %Geometry Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = public)
        type
        cartd
        cart_pos_gp
        shape
        dvolu
        djacob
        interpolation
        quadrature
        nfields=1;
    end
    methods (Access = {?Physical_Problem,?Element_DiffReact}) % !! Element_DiffReact -> Chapusilla !!
        function obj = Geometry(mesh,interp_order,order)
            if nargin==2
                order=interp_order;
            end    
            obj.type=mesh.geometryType;
            obj.interpolation=Interpolation(mesh,interp_order);            
            obj.computeGeometry(order);
        end
        function computeGeometry(obj,order)
            obj.createQuadrature(order);
            ndime=obj.interpolation.isoparametric.ndime;
            nnode=obj.interpolation.isoparametric.nnode;
            ngaus=obj.quadrature.ngaus;
            nelem=obj.interpolation.nelem;
            gp_position = zeros(ndime,obj.quadrature.ngaus,nelem);
            
            obj.shape=zeros(nnode,ngaus);
            obj.cartd=zeros(ndime,nnode,nelem,ngaus);
            for igauss=1:ngaus
                pos_gp = num2cell(obj.quadrature.posgp(:,igauss));
                obj.shape(:,igauss) = cell2mat(obj.interpolation.isoparametric.shape(pos_gp{:}));
                
                for inode = 1:nnode
                    for idime=1:ndime
                        x = obj.shape(inode,igauss).*squeeze(obj.interpolation.xpoints(obj.interpolation.T(:,inode),idime));
                        %                         y = shape(2,igauss).*squeeze(interpolation.xpoints(interpolation.T(inode),2));
                        gp_position(idime,igauss,:) = squeeze(gp_position(idime,igauss,:))+x;
                    end
                end
            end
            
            %             end
            
            obj.cart_pos_gp = gp_position;
            for i = 1:ndime
                a = obj.interpolation.xpoints(:,i);
                elcoord(:,i,:) = a(permute(obj.interpolation.T',[1,3,2]));
            end
            % Gauss loop
            for igauss = 1:ngaus
                pos_gp = num2cell(obj.quadrature.posgp(1:ndime,igauss));
                % shape(:,igauss) = cell2mat(interpolation.isoparametric.shape(pos_gp{:}));
                deriv(:,:,igauss) = cell2mat (obj.interpolation.isoparametric.deriv(pos_gp{:}));
                for i = 1:ndime
                    % Jacobian
                    deriv_perm = permute(deriv(i,:,igauss),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,ndime,nelem);
                    jacobian(i,:,:) = sum(deriv_perm_large.*elcoord,1);
                end
                
                [invJ,detJ]=obj.inverseElementalMatrix(ndime,jacobian);
                
                for i = 1:ndime
                    % Cartesian Derivatives
                    deriv_perm = permute(invJ(i,:,:),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,obj.interpolation.isoparametric.nnode,1) .*repmat(deriv(:,:,igauss),1,1,nelem);
                    obj.cartd(i,:,:,igauss) = sum(deriv_perm_large);
                end
                %                 w= num2cell(weigp(igauss));
                %                 obj.weigp(igauss) =  cell2mat(w{igauss});
                obj.dvolu(:,igauss) = obj.quadrature.weigp(igauss)*detJ;
                obj.djacob(:,igauss)= detJ;
            end
            % dvolu         
            
        end
        function createQuadrature(obj,order)            
            switch obj.type             
                case 'TRIANGLE'
                    switch order
                        case 'LINEAR'
                            obj.quadrature = Quadrature_Triangle(order);
                        case 'QUADRATIC'
                            obj.quadrature = Quadrature_Triangle(order);
                        otherwise
                            error('Invalid order for quadrature TRIANGLE.');
                    end
                case 'QUAD'
                    switch order
                        case 'CONSTANT'
                            
                        case 'LINEAR'
                            obj.quadrature = Quadrature_Quadrilateral(order);
                        case 'QUADRATIC'
                            obj.quadrature = Quadrature_Quadrilateral(order);
                        otherwise
                            error('Invalid order for quadrature QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                    obj.quadrature = Quadrature_Tetrahedra(order);
                case 'HEXAHEDRA'
                    obj.quadrature = Quadrature_Hexahedra(order);
                otherwise
                    error('Invalid mesh type.')
            end
        end
    end
    methods (Static)
        function [inverse,determinant]=inverseElementalMatrix(ndime,A)
            switch ndime
                case 2
                    [inverse,determinant] = multinverse2x2(A);
                case 3
                    [inverse,determinant] = multinverse3x3(A);
            end
        end
        
    end
    
end