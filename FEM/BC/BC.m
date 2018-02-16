classdef BC < handle
    %BC Summary of this class goes here
    %   Detailed explanation goes here
    
    

    properties (GetAccess = public, SetAccess = public)
        fixnodes
        fixnodes_perimeter
        neunodes
        iN
        iD
    end
    
    methods (Access = public)
        % Constructor
        function obj = BC(nunkn,filename)
            [obj.fixnodes,obj.fixnodes_perimeter,obj.neunodes] = Preprocess.getBC(filename);
            obj.computeiDiN(nunkn);
        end
        
        function  obj = create(nunkn,filename)
            
            
        end
    end
    
    methods (Access = protected)
        function obj = computeiDiN(obj,nunkn)
            obj.iD = obj.comupute_global_id_nodes_from_local(obj.fixnodes,nunkn);
            obj.iN = obj.comupute_global_id_nodes_from_local(obj.neunodes,nunkn);
        end        
    end
    
    methods (Static)
        function iG = comupute_global_id_nodes_from_local(iL,nunkn)
            
            if (~isempty(iL))
                iG = zeros(size(iL,1),1);
                for i = 1:length(iL(:,1))
                    iG(i,1)= nunkn*(iL(i,1) - 1) + iL(i,2);
                end
            end
        end
    end
end
    
