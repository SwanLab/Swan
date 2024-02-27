classdef UnfittedFunction < handle

    properties (Access = public)
        ndimf
        %nDofs
        %order
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
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.unfittedMesh = cParams.uMesh;
            obj.fun          = cParams.fun;
            obj.ndimf        = cParams.fun.ndimf;
        end

        function computeUnfittedMeshFunction(obj)
            uMeshFun  = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(obj.fun);
            %obj.nDofs = uMeshFun.backgroundFunction.nDofs;
            %obj.order = uMeshFun.backgroundFunction.order;
            obj.innerMeshFunction    = uMeshFun.innerMeshFunction;
            obj.innerCutMeshFunction = uMeshFun.innerCutMeshFunction;
        end
    end
end