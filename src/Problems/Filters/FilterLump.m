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
            xFun.setFValues(full(xProj));
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
            f   = @(v,u) DP(v,u);
            lhs = IntegrateLHS(f,obj.trial,obj.trial,obj.mesh,'Domain',2);
            obj.LHS = obj.lumpMatrix(lhs);
        end

        function RHS = computeRHS(obj,fun,quadType)
            f   = @(v) DP(fun,v);
            RHS = IntegrateRHS(f,obj.trial,obj.mesh,'Domain',quadType);
        end

    end

    methods (Access = private, Static)
        function Al = lumpMatrix(A)
            I  = ones(size(A,2),1);
            Al = A*I;
        end
    end

end
