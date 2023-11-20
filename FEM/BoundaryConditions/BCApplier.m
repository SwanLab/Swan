classdef BCApplier < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = BCApplier(cParams)
            obj.init(cParams)
        end
        
    end

    methods (Static, Access = public)

        function HSr = reduceLHS(HS, bcs)
            cstr = BCApplier.getConstrainedNodes(bcs);
            ndofs = 1:length(HS);
            free = setdiff(ndofs, cstr);

            if (size(HS,1) == size(HS,2)) % Matrix
                HSr = HS(free, free);
            else %Vector
                HSr = HS(free);
            end
        end
    end
    
    methods (Static, Access = private)
        
        function dofs = getConstrainedNodes(bcs)
            dofs = [];
            for iBc = 1:length(bcs)
                bc = bcs{iBc};
                if ~strcmp(bc.type, 'Neumann')
                    a = 1;
                    dofs = [dofs; bc.getDofs];
                end
            end
        end
        
    end

    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end