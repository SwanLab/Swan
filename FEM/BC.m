classdef BC
    %BC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Physical_Problem,?Element,?Filter}, SetAccess = private)
        fixnodes
        fixnodes_perimeter
        neunodes
        iN
        iD
    end
    
    methods (Access = ?Physical_Problem)
        % Constructor
        function obj = BC(nunkn,filename)
            if nargin ~= 0 
                [obj.fixnodes,obj.fixnodes_perimeter,obj.neunodes] = Preprocess.getBC(filename);
                
                for i = 1:length(obj.fixnodes(:,1))
                    obj.iD(i) = obj.fixnodes(i,1)*nunkn - nunkn + obj.fixnodes(i,2);
                end
                
                for i = 1:length(obj.neunodes(:,1))
                    obj.iN(i)= obj.neunodes(i,1)*nunkn - nunkn + obj.neunodes(i,2);
                end
                
            end
        end
    end
    
end

