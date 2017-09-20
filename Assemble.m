classdef Assemble
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % test
    end
    
    methods (Static)
        function [LHS,RHS] = Compute(element,nnode,nunkn,dof)
            
            % Compute LHS
            
            LHS = sparse(dof.ndof,dof.ndof);
            
            for i = 1:nnode*nunkn
                for j = 1:nnode*nunkn
                    a = squeeze(element.LHS(i,j,:));
                    LHS = LHS + sparse(dof.idx(i,:),dof.idx(j,:),a,dof.ndof,dof.ndof);
                end
            end
            
            % Compute RHS
            
            RHS = zeros(1,dof.ndof);
            for i = 1:length(dof.idx(:,1)) % nnode*nunkn
               b = squeeze(element.RHS(i,1,:));
               ind = dof.idx(i,:);
               RHS(ind) = b; 
            end
            
        end
    end
    
end

