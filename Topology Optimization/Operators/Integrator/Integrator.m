classdef Integrator < handle

    methods (Access = public, Static)
        function obj = create(type,mesh,order)
            s.type     = type;
            s.mesh     = mesh;
            s.quadType = order;
            f          = IntegratorFactory();
            obj        = f.create(s);
        end
    end

end