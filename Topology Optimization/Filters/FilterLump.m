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
            nVar = numel(x);
            xFun = cell(nVar,1);
            for ivar = 1:nVar
                xVar         = x{ivar};
                xFun{ivar}   = LagrangianFunction.create(obj.mesh, xVar.ndimf, obj.trial.order);
                lhs          = obj.LHS;
                rhs          = obj.computeRHS(xVar, quadType);
                xProj        = rhs./lhs;
                xFun{ivar}.fValues = xProj;
            end
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial         = cParams.trial.copy();
            obj.mesh          = cParams.mesh;
        end

        function computeLHS(obj)
            s.mesh            = obj.mesh;
            s.test            = obj.trial;
            s.trial           = obj.trial;
            s.quadratureOrder = 'QUADRATIC';
            s.type            = 'MassMatrix';
            int               = LHSintegrator.create(s);
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
            int        = RHSintegrator.create(s);
            test       = LagrangianFunction.create(obj.mesh,fun.ndimf,obj.trial.order);
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