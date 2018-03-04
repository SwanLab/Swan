classdef Geometry
    %Geometry Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        cartd
        cart_pos_gp
        shape
        dvolu
        weigp
        ndime
        nnode
        ngaus
        djacob
                
    end
    properties (GetAccess = ?Postprocess, SetAccess = private)
        posgp
    end
    
    methods (Access = ?Physical_Problem)
        function obj = Geometry(interpolation,quadrature,nelem)
            
            geometryObject = interpolation.isoparametric;
            obj.ndime = interpolation.isoparametric.ndime;
            obj.posgp = quadrature.posgp;
            weigp = quadrature.weigp;
            obj.ngaus = quadrature.ngaus;
%             obj.shape = interpolation.isoparametric.shape;
            obj.nnode = interpolation.isoparametric.nnode;
            ndime = interpolation.isoparametric.ndime;
            
            gp_position = zeros(obj.ndime,obj.ngaus,nelem);    
%             for ielem=1:nelem
%                 T_elem= interpolation.T(ielem,:);
                for igauss=1:obj.ngaus
                    pos_gp = num2cell(obj.posgp(1:obj.ndime,igauss));
                    shape(:,igauss) = cell2mat(interpolation.isoparametric.shape(pos_gp{:}));
                    
                    for inode = 1:obj.nnode
                        for idime=1:ndime
                        x = shape(inode,igauss).*squeeze(interpolation.xpoints(interpolation.T(:,inode),idime));
%                         y = shape(2,igauss).*squeeze(interpolation.xpoints(interpolation.T(inode),2));
                        gp_position(idime,igauss,:) = squeeze(gp_position(idime,igauss,:))+x;  
                        end
                    end
                end
                
%             end
            
                obj.cart_pos_gp = gp_position;

            for i = 1:obj.ndime
                a = interpolation.xpoints(:,i);
                elcoord(:,i,:) = a(permute(interpolation.T',[1,3,2]));
            end
            % Gauss loop
            for igauss = 1:obj.ngaus
                 pos_gp = num2cell(obj.posgp(1:obj.ndime,igauss));
                 obj.shape(:,igauss) = cell2mat(interpolation.isoparametric.shape(pos_gp{:}));
                 deriv(:,:,igauss) = cell2mat (interpolation.isoparametric.deriv(pos_gp{:}));
                for i = 1:obj.ndime
                    % Jacobian
                    deriv_perm = permute(deriv(i,:,igauss),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,obj.ndime,nelem);
                    jacobian(i,:,:) = sum(deriv_perm_large.*elcoord,1);
                end                
                % !! SWITCH EXECPCIÓ? !!
                switch obj.ndime
                    case 2
                        [invJ,detJ] = multinverse2x2(jacobian);
                    case 3
                        [invJ,detJ] = multinverse3x3(jacobian);
                end
                
                for i = 1:obj.ndime
                    % Cartesian Derivatives
                    deriv_perm = permute(invJ(i,:,:),[2,1,3]);
                    deriv_perm_large = repmat(deriv_perm,1,geometryObject.nnode,1) .*repmat(deriv(:,:,igauss),1,1,nelem);
                    obj.cartd(i,:,:,igauss) = sum(deriv_perm_large);
                end
                w= num2cell(weigp(igauss));
                obj.weigp(igauss) =  cell2mat(w{igauss});
                obj.dvolu(:,igauss) = obj.weigp(igauss)*detJ;
                obj.djacob(:,igauss)= detJ;
            end
            % dvolu
            
            
        end
    end
    
end