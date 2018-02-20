classdef Element_Elastic_2D<Element_Elastic
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function variables = computeVars(obj,uL)
            variables = obj.computeDispStressStrain(uL);
            variables.strain = obj.computeEz(variables.strain,obj.nstre,obj.nelem,obj.material);
            variables = obj.permuteStressStrain(variables);
        end
        
        
        function [B] = computeB(obj,nunkn,nelem,nnode,cartd)
            B = zeros(3,nnode*nunkn,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)  = cartd(1,i,:);
                B(2,j+1,:)= cartd(2,i,:);
                B(3,j,:)  = cartd(2,i,:);
                B(3,j+1,:)= cartd(1,i,:);
            end
        end
    end
    
end

