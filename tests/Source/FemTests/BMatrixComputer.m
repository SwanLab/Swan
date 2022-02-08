classdef BMatrixComputer < handle
    % Ideally, this class should disappear soon.
    
    properties (Access = private)
        dim
        geometry
    end
    
    methods (Access = public)

        function obj = BMatrixComputer(cParams)
            obj.init(cParams)
        end

        function B = compute(obj, igaus)
            B = obj.computeB(igaus);
        end

    end

    methods (Access = private)
        
        function init(obj, cParams)
            obj.dim = cParams.dim;
            obj.geometry = cParams.geometry;
        end

       function B = computeB(obj,igaus)
           ndim = obj.dim.ndim;
           switch ndim
               case 2
                   B = obj.computeB2D(igaus);
               case 3
                   B = obj.computeB3D(igaus);
           end
       end

       function B = computeB2D(obj,igaus)
            nstre = obj.dim.nstre;
            nnode = obj.dim.nnode;
            nelem = obj.dim.nelem;
            nunkn = obj.dim.nunkn; 
            ndofPerElement = obj.dim.ndofPerElement;
            B = zeros(nstre,ndofPerElement,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)  = obj.geometry.cartd(1,i,:,igaus);
                B(2,j+1,:)= obj.geometry.cartd(2,i,:,igaus);
                B(3,j,:)  = obj.geometry.cartd(2,i,:,igaus);
                B(3,j+1,:)= obj.geometry.cartd(1,i,:,igaus);
            end
       end

       function [B] = computeB3D(obj,igaus)
           d = obj.dim;
            B = zeros(d.nstre,d.ndofPerElement,d.nelem);
            for inode=1:d.nnode
                j = d.nunkn*(inode-1)+1;
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

