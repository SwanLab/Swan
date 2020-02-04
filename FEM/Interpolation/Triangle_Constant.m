classdef Triangle_Constant < Interpolation

    methods (Access = public)
        function obj = Triangle_Constant(mesh)
            obj = obj@Interpolation(mesh);
            obj.type = 'TRIANGLE';
            obj.ndime = 2;
            obj.nnode = 1;
            obj.pos_nodes = [1/3 1/3];
            obj.shape = @(s,t) {1};
            obj.deriv = @(s,t) {0;0};
        end
    end
end
