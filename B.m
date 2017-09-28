classdef B<Matrix_Local
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = B(nstre,nunkn,nelem,nnode,cartd,ndime)
            obj.value = zeros(nstre,nnode*nunkn,nelem);
            switch ndime
                case 2
                    for i = 1:nnode
                        j = nunkn*(i-1)+1;
                        obj.value(1,j,:)  = cartd(1,i,:);
                        obj.value(2,j+1,:)= cartd(2,i,:);
                        obj.value(3,j,:)  = cartd(2,i,:);
                        obj.value(3,j+1,:)= cartd(1,i,:);
                    end
                case 3
                    for inode=1:nnode
                        j = nunkn*(inode-1)+1;
                        % associated to normal strains
                        obj.value(1,j,:)=cartd(1,inode,:);
                        obj.value(2,j+1,:)=cartd(2,inode,:);
                        obj.value(3,j+2,:)=cartd(3,inode,:);
                        % associated to shear strain, gamma12
                        obj.value(4,j,:)=cartd(2,inode,:);
                        obj.value(4,j+1,:)=cartd(1,inode,:);
                        % associated to shear strain, gamma13
                        obj.value(5,j,:)=cartd(3,inode,:);
                        obj.value(5,j+2,:)=cartd(1,inode,:);
                        % associated to shear strain, gamma23
                        obj.value(6,j+1,:)=cartd(3,inode,:);
                        obj.value(6,j+2,:)=cartd(2,inode,:);
            end
            end
        end
    end
    
end

