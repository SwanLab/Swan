classdef Element_Elastic<Element
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nstre = 3;
    end
    
    
    methods
        function obj = computeLHS(obj,nunkn,nelem,geometry)
            Kel = Ke(obj.nstre,nunkn,nelem,geometry);
            obj.LHS = Kel.value;
        end
        
        function obj = computeRHS(obj,nunkn,nelem,nnode,bc,idx)
            Fext = zeros(nnode*nunkn,1,nelem);
            for i = 1:length(bc.iN)
                for j = 1:nelem
                    ind = find(idx(:,j) == bc.iN(i));
                    if ~isempty(ind)
                        Fext(ind,:,j) = bc.Fpointload(i,3);
                    end
%                     clear ind
                     ind = [];
                end
            end
            obj.RHS = Fext;
        end
        
    end
end
