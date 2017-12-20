classdef BC
    %BC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Physical_Problem,?Element}, SetAccess = protected)
        fixnodes
        neunodes
        iN
        iD
    end
    
    methods (Access = public)
        % Constructor
        function obj = BC(nunkn,filename)
            [obj.fixnodes,obj.neunodes] = Preprocess.getBC(filename);
            
            if (~isempty(obj.fixnodes))
                for i = 1:length(obj.fixnodes(:,1))
                    obj.iD(i) = obj.fixnodes(i,1)*nunkn - nunkn + obj.fixnodes(i,2);
                end
            end
            if (~isempty(obj.neunodes))
                for i = 1:length(obj.neunodes(:,1))
                    obj.iN(i)= obj.neunodes(i,1)*nunkn - nunkn + obj.neunodes(i,2);
                end
            end
        end
    end
    
end

