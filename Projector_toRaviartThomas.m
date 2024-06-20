classdef Projector_toRaviartThomas < Projector
    
    properties (Access = private)
        order
    end
    
    methods (Access = public)

        function obj = Projector_toRaviartThomas(cParams)
            obj.init(cParams);
            obj.order = cParams.projectorType;
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS(x);
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = reshape(xProj,[x.ndimf/2,numel(xProj)/x.ndimf*2])'; % no
            s.order = obj.order;
            xFun = RaviartThomasFunction(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj,fun)
            s.mesh  = obj.mesh;
            s.test  = RaviartThomasFunction.create(obj.mesh, fun.ndimf/2, obj.order);
            s.trial = RaviartThomasFunction.create(obj.mesh, fun.ndimf/2, obj.order);
            % s.quadratureOrder = 'QUADRATIC'; % no
            s.type  = 'MassMatrixRT';
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
            test       = RaviartThomasFunction.create(obj.mesh,fun.ndimf/2,obj.order);
            RHS        = int.compute(fun,test);
        end

        function ord = createRHSQuadrature(obj, fun)
            ord = 2;
        end

    end

end