classdef UnfittedBoundaryFunction < handle

    properties (Access = public)
        ndimf
        unfittedMesh
        boundaryCutMeshFunction
        unfittedBoundaryMeshFunction
    end

    properties (Access = private)
        fun
    end

    methods (Access = public)
        function obj = UnfittedBoundaryFunction(cParams)
            obj.init(cParams);
            obj.computeUnfittedMeshFunction();
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
            obj.boundaryCutMeshFunction = uMeshFun.boundaryCutMeshFunction;
            obj.unfittedBoundaryMeshFunction = uMeshFun.unfittedBoundaryMeshFunction;
        end
    end
end