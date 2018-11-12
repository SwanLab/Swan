classdef Line_Linear < handle
    properties (GetAccess = public, SetAccess = private)
        ndime
        nnode
        order
        npnod
        nelem
        type
        pos_nodes
        shape
        deriv
        dvolu
    end
    
    methods (Access = public)
        function obj = Line_Linear
            obj.type = 'LINE';
            obj.order = 'LINEAR';
            obj.ndime = 1;
            obj.nnode = 2;
            obj.pos_nodes = [-1; 1];
            obj.dvolu = 2;
        end
        
        function computeShapeDeriv(obj,posgp)
            obj.shape = [];
            obj.deriv = [];
            s = posgp;
            
            obj.shape = [ones(length(posgp),1)-s,s+1]/2;
            obj.deriv = repmat([-1.0,1.0],1,1,length(posgp));
        end
    end
end
