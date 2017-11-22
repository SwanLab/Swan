classdef B_thermal<B
    %B2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Access = ?Element_Thermal)
        function [obj,B] = computeB(obj,nunkn,nelem,nnode,cartd)
            B = zeros(2,nnode*nunkn,nelem);
            for inode=1:nnode
                j = nunkn*(inode-1)+1;
                B(1,j,:)=cartd(1,inode,:);
                B(2,j,:)=cartd(2,inode,:);
            end
            obj.value = [obj.value {B}];
        end
    end
    
end
