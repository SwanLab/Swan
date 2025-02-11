classdef ConnecRenumbering < handle
       
    properties (Access = private)
        oldNodes
        newNodes
    end
    
    methods (Access = public)
        
        function obj = ConnecRenumbering(cParams)
            obj.init(cParams)
        end

        function newConnec = renumber(obj,connec)
            nodes        = obj.oldNodes;
            M            = zeros(max(nodes(:)),1);
            M(nodes)     = obj.newNodes;
            newConnec    = zeros(size(connec));
            newConnec(:) = M(connec);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.oldNodes = cParams.oldNodes;
            obj.newNodes = cParams.newNodes;
        end
        
    end
    
end