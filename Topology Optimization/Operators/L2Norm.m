classdef L2Norm < handle

    properties (Access = private)
        integrator
    end

    methods (Access = public)
        function obj = L2Norm(mesh)
            obj.createIntegrator(mesh);
        end

        function sp = compute(obj,f)
            sp = obj.integrator.compute(f,f);
        end
    end

    methods (Access = private)
        function createIntegrator(obj,m)
            q.mesh         = m;
            q.quadType     = 'QUADRATIC';
            q.type         = 'ScalarProduct';
            int            = Integrator.create(q);
            obj.integrator = int;
        end
    end
end