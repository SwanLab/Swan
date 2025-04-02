classdef UnfittedFunction < handle

    properties (Access = public)
        ndimf
        unfittedMesh
        innerMeshFunction
        innerCutMeshFunction
    end

    properties (Access = private)
        fun
    end

    methods (Access = public)
        function obj = UnfittedFunction(cParams)
            obj.init(cParams);
            obj.computeUnfittedMeshFunction();
        end

        function f = copy(obj)
            s.uMesh = obj.unfittedMesh;
            s.fun   = obj.fun;
            s.ndimf = obj.ndimf;
            f       = UnfittedFunction(s);
        end

        function res = times(obj1,obj2)
            res     = copy(obj1);
            res.fun = res.fun.*obj2;
            res.computeUnfittedMeshFunction();
        end

        function fun = project(obj,target)
            s.mesh          = obj.unfittedMesh.backgroundMesh;
            s.projectorType = target;
            proj            = Projector.create(s);
            fun             = proj.project(obj);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.unfittedMesh = cParams.uMesh;
            obj.fun          = cParams.fun;
            obj.ndimf        = cParams.fun.ndimf;
        end

        function computeUnfittedMeshFunction(obj)
            uMeshFun = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(obj.fun);
            obj.innerMeshFunction    = uMeshFun.innerMeshFunction;
            obj.innerCutMeshFunction = uMeshFun.innerCutMeshFunction;
        end
    end
end