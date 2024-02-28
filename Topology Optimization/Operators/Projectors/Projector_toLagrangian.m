classdef Projector_toLagrangian < Projector
    
    properties (Access = private)
        order
    end
    
    methods (Access = public)

        function obj = Projector_toLagrangian(cParams)
            obj.init(cParams);
            obj.order = cParams.projectorType;
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            s.order = obj.order;
            xFun = LagrangianFunction(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj)
            s.mesh  = obj.mesh;
            s.test  = LagrangianFunction.create(obj.mesh, 1, obj.order);
            s.trial = LagrangianFunction.create(obj.mesh, 1, obj.order);
            s.quadratureOrder = 'ORDER10'; % no
            s.type  = 'MassMatrix';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
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
            test       = LagrangianFunction.create(obj.mesh,1,obj.order);
            RHS        = int.compute(fun,test);
        end

        function ord = createRHSQuadrature(obj, fun)
            if isa(fun, 'FGaussDiscontinuousFunction')
                ord = fun.getQuadratureOrder();
            else
                ord = obj.determineQuadratureOrder(fun);
                ord = 'ORDER10'; % no
                ord = 'QUADRATIC';
            end
        end

    end

end