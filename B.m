classdef B<Matrix_Local
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = B(nstre,nunkn,nelem,nnode,cartd)
            obj.value = zeros(nstre,nnode*nunkn,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                obj.value(1,j,:)  = cartd(1,i,:);
                obj.value(2,j+1,:)= cartd(2,i,:);
                obj.value(3,j,:)  = cartd(2,i,:);
                obj.value(3,j+1,:)= cartd(1,i,:);
            end
        end
    end
    
end

