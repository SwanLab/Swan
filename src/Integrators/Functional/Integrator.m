classdef Integrator < handle

    methods (Access = public, Static)
        function int = compute(fun,mesh,order)
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    type = 'Unfitted';
                    m    = fun.unfittedMesh;
                otherwise
                    type = 'Function';
                    m    = mesh;
            end
            integrator = Integrator.create(type,m,order);
            int        = integrator.compute(fun);
        end

        function obj = create(type,mesh,order)
            s.type     = type;
            s.mesh     = mesh;
            s.quadType = order;
            f          = IntegratorFactory();
            obj        = f.create(s);
        end
    end

end