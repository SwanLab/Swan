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
            LHS = obj.computeLHS(x);
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = reshape(xProj,[x.ndimf,numel(xProj)/x.ndimf])';
            s.order = obj.order;
            xFun = LagrangianFunction(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj,fun)
            s.mesh  = obj.mesh;
            s.test  = LagrangianFunction.create(obj.mesh, fun.ndimf, obj.order);
            s.trial = LagrangianFunction.create(obj.mesh, fun.ndimf, obj.order);
            s.quadratureOrder = 2; % no
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
            test       = LagrangianFunction.create(obj.mesh,fun.ndimf,obj.order);
            RHS        = int.compute(fun,test);
        end

        function ord = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
        end

    end

end