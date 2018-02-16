classdef DOF_micro < DOF
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vP
        vQ
    end
    
    methods
        
        % Constructor
        function obj = DOF_micro()
            
        end
        
       function vC = compute_constrained_nodes(obj)
            vC = [obj.vQ;obj.vD];
       end
        
    end
    
    methods (Static)
        function [vp,vq] = compute_periodic_nodes(nunkn,bc)
            nlib = size(bc.pnodes(1,:),2);
            vp = zeros(nlib*nunkn,1);
            vq = zeros(nlib*nunkn,1);
            for iunkn = 1:nunkn
                index_glib = nlib*(iunkn - 1) + [1:nlib];
                vp(index_glib,1) = (bc.pnodes(1,:)-1)*nunkn + iunkn;
                vq(index_glib,1) = (bc.pnodes(2,:)-1)*nunkn + iunkn;
            end
        end
            
        
        
    end
    
end

