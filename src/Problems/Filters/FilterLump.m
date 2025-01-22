classdef FilterLump < handle

    properties (Access = private)
        trial
        mesh
    end

    properties (Access = private)
        LHS
    end

    methods (Access = public)

        function obj = FilterLump(cParams)
            obj.init(cParams);
            obj.computeLHS();
        end

        function xFun = compute(obj, x, quadType)
            s.feFunType  = class(obj.trial);
            s.mesh       = obj.mesh;
            s.ndimf      = 1;
            xFun         = LagrangianFunction.create(s.mesh, s.ndimf, obj.trial.order);
            lhs          = obj.LHS;
            rhs          = obj.computeRHS(x, quadType);
            xProj        = rhs./lhs;
            xFun.setFValues(xProj);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            cParams.feFunType = class(cParams.trial);
            cParams.ndimf     = 1;
            obj.trial         = LagrangianFunction.create(cParams.mesh, 1, cParams.trial.order);
            obj.mesh          = cParams.mesh;
        end

        function computeLHS(obj)
            s.mesh            = obj.mesh;
            s.test            = obj.trial;
            s.trial           = obj.trial;
            s.quadratureOrder = 2;
            s.type            = 'MassMatrix';
            int               = LHSIntegrator.create(s);
            lhs               = int.compute();
            obj.LHS           = obj.lumpMatrix(lhs);
        end

        function rhs = computeRHS(obj,fun,quadType)
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    s.mesh = fun.unfittedMesh;
                    s.type = 'Unfitted';
                otherwise
                    s.mesh = obj.mesh;
                    s.type = 'ShapeFunction';
            end
            s.quadType = quadType;
            int        = RHSIntegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(fun,test);
        end

    end

    methods (Access = private, Static)
        function Al = lumpMatrix(A)
            I  = ones(size(A,2),1);
            Al = A*I;
        end
    end

end