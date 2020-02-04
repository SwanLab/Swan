classdef Line_Linear < Interpolation
    
    methods (Access = public)
        function obj = Line_Linear(cParams)
            obj.init(cParams);
            obj.type = obj.mesh.geometryType;
            obj.ndime = 1;
            obj.nnode = 2;
            obj.pos_nodes = [-1; 1];
            obj.isoDv = 2;
            obj.computeCoordAndConnec();
        end
        
        function computeShapeDeriv(obj,posgp)
            ngaus = size(posgp,2);
            s = posgp(1,:);
            I = ones(1,ngaus);
            obj.shape = [I-s;s+1]/2;
            obj.deriv = repmat([-0.5,0.5],1,1,ngaus);
        end
        
    end
    
end
