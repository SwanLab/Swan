classdef BC_Micro < BC
    %BC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Physical_Problem,?Element}, SetAccess = private)
        pnodes
    end
    
    methods (Access = ?Physical_Problem)
        % Constructor
        function obj = BC_Micro(nunkn,filename,coords)
            obj@BC(nunkn,filename);
            [obj.fixnodes,obj.pnodes] = Preprocess.getPeriodicBC(coords);      
        end
    end
end


