classdef Projector_toNedelec < Projector
    
    properties (Access = private)
        order
    end
    
    methods (Access = public)

        function obj = Projector_toNedelec(cParams)
            obj.init(cParams);
            obj.order = cParams.projectorType;
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS(x);
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = reshape(xProj,[1,numel(xProj)])'; % no
            s.order = obj.order;
            xFun = NedelecFunction(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj,fun)
            s.mesh  = obj.mesh;
            s.test  = NedelecFunction.create(obj.mesh, 1, obj.order);
            s.trial = NedelecFunction.create(obj.mesh, 1, obj.order);
            % s.quadratureOrder = 'QUADRATIC'; % no
            s.type  = 'MassMatrixVect';
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
                    s.type = 'ShapeFunctionN';
            end
            s.quadType = ord;
            int        = RHSintegrator.create(s);
            test       = NedelecFunction.create(obj.mesh,1,obj.order);
            RHS        = int.compute(fun,test);
        end

        function ord = createRHSQuadrature(obj, fun)
            ord = 2;
        end

    end

end