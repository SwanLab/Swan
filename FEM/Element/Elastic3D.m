classdef Elastic3D < ElasticDim

    
    properties (Access = protected)
        nstre = 6;
    end
    
    methods (Access = public)
                
        function [B] = computeB(obj,igaus)
            B = zeros(obj.nstre,obj.nnode*obj.dof.nunkn,obj.nelem);
            for inode=1:obj.nnode
                j = obj.dof.nunkn*(inode-1)+1;
                % associated to normal strains
                B(1,j,:) = obj.geometry.cartd(1,inode,:,igaus);
                B(2,j+1,:) = obj.geometry.cartd(2,inode,:,igaus);
                B(3,j+2,:) = obj.geometry.cartd(3,inode,:,igaus);
                % associated to shear strain, gamma12
                B(4,j,:) = obj.geometry.cartd(2,inode,:,igaus);
                B(4,j+1,:) = obj.geometry.cartd(1,inode,:,igaus);
                % associated to shear strain, gamma13
                B(5,j,:) = obj.geometry.cartd(3,inode,:,igaus);
                B(5,j+2,:) = obj.geometry.cartd(1,inode,:,igaus);
                % associated to shear strain, gamma23
                B(6,j+1,:) = obj.geometry.cartd(3,inode,:,igaus);
                B(6,j+2,:) = obj.geometry.cartd(2,inode,:,igaus);
            end
        end
                
    
        
    end
end



