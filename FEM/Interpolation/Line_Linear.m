classdef Line_Linear < Interpolation

    methods (Access = public)
        function obj = Line_Linear(mesh)
            obj = obj@Interpolation(mesh);
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
            obj.deriv = repmat([-0.5,0.5],1,1,length(posgp));
        end
    end
end
