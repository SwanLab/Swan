classdef Integrator < handle

    methods (Access = public, Static)
        function obj = create(s)
            f   = IntegratorFactory();
            obj = f.create(s);
        end
    end

end