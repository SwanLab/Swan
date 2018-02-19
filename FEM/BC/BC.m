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
    end
    
    methods (Access = public)
        function obj = computeiDiN(obj,nunkn)
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
    
