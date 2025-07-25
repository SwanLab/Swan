classdef UnfittedFunction < BaseFunction

    properties (Access = public)
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
            obj.mesh = obj.unfittedMesh.backgroundMesh;
        end

        function res = times(obj1,obj2)
            res     = copy(obj1);
            res.fun = res.fun.*obj2;
            res.computeUnfittedMeshFunction();
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

    methods (Access = protected)

        function evaluateNew(obj,xV)

        end

    end
end