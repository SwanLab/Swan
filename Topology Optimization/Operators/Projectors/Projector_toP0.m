classdef Projector_toP0 < Projector

    properties (Access = private)
        quadrature
        M
    end

    methods (Access = public)

        function obj = Projector_toP0(cParams)
            obj.init(cParams);
            obj.computeQuadrature();
            obj.createMassMatrix();
        end

        function xFun = project(obj, x)
            RHS = obj.createRHS(x);
            s.fValues = obj.M\RHS;
            s.mesh    = obj.mesh;
            s.order   = 'P0';
            xFun = LagrangianFunction(s);
        end

    end

    methods (Access = private)

        function createMassMatrix(obj)
            quad = Quadrature.create(obj.mesh,1);
            dv = obj.mesh.computeDvolume(quad);
            a = sum(dv(1,:),1);
            obj.M = spdiags(a',0,length(a),length(a));
         %   obj.M = spdiags(sum(dv(1,:),1),0);
        end

        function rhs = createRHS(obj, fun)
            ord = obj.createRHSQuadrature(fun);
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    s.mesh = fun.unfittedMesh;
                    s.type = 'Unfitted';
                otherwise
                    s.mesh = obj.mesh;
                    s.type = 'ShapeFunction';
            end
            s.quadType = ord;
            int        = RHSintegrator.create(s);
            test       = LagrangianFunction.create(obj.mesh,fun.ndimf,'P0');
            rhs        = int.compute(fun,test);
        end

        function ord = createRHSQuadrature(obj, fun)
            if isa(fun, 'FGaussDiscontinuousFunction')
                ord = fun.getQuadratureOrder();
            else
                ord = obj.determineQuadratureOrder(fun);
                ord = 'QUADRATIC'; % no
            end
        end

        function computeQuadrature(obj)
            quad = Quadrature.create(obj.mesh,2);
            obj.quadrature = quad;
        end
    end

end