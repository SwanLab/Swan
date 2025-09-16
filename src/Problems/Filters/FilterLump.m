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
            s.ndimf      = x.ndimf;
            xFun         = LagrangianFunction.create(s.mesh, s.ndimf, obj.trial.order);
            lhs          = obj.LHS;
            rhs          = obj.computeRHS(x, quadType);
            xProj        = rhs./lhs;
            xProj        = reshape(xProj',obj.trial.ndimf,[])';
            xFun.setFValues(xProj);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            cParams.feFunType = class(cParams.trial);
            cParams.ndimf     = cParams.trial.ndimf;
            obj.trial         = LagrangianFunction.create(cParams.mesh, cParams.ndimf, cParams.trial.order);
            obj.mesh          = cParams.mesh;
        end

        function computeLHS(obj)
            s.mesh  = obj.trial.mesh;
            s.quadratureOrder = 2; % no
            lhs = LHSIntegrator(s);
            for i = 1:obj.trial.nDofsElem
                v = Test(obj.trial,i);
                for j = 1:obj.trial.nDofsElem
                    u = Test(obj.trial,j);
                    f{i,j} = DP(v,u);
                end
            end
            lhs = lhs.compute(f,obj.trial,obj.trial);
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
