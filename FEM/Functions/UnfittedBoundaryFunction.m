classdef UnfittedBoundaryFunction < handle

    properties (Access = public)
        unfittedMesh
    end

    properties (Access = private)
        ndimf
        fun
    end

    properties (Access = private)
        unfittedMeshFunction
    end

    methods (Access = public)

        function obj = UnfittedBoundaryFunction(cParams)
            obj.init(cParams);
            obj.computeUnfittedMeshFunction();
        end

        function f = obtainFunctionAtExternalBoundary(obj,iBoundary)
            uMesh      = obj.unfittedMesh;
            meshes     = uMesh.unfittedBoundaryMesh.getActiveMesh();
            uMeshFun   = obj.unfittedMeshFunction;
            uBMeshFuns = uMeshFun.unfittedBoundaryMeshFunction;
            s.uMesh    = meshes{iBoundary};
            s.fun      = uBMeshFuns.activeFuns{iBoundary}.backgroundFunction;
            f          = UnfittedFunction(s);
        end

        function fxV = evaluateCutElements(obj,q)
            uMeshFun = obj.unfittedMeshFunction;
            fNew     = uMeshFun.boundaryCutMeshFunction;
            xV       = q.posgp;
            fxV      = fNew.evaluate(xV);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.unfittedMesh   = cParams.uMesh;
            obj.fun            = cParams.fun;
            obj.ndimf          = cParams.fun.ndimf;
        end

        function computeUnfittedMeshFunction(obj)
            fType = obj.fun.fType;
            switch fType
                case 'L2' % Provisional
                    f = obj.fun.project('P1');
                otherwise
                    f = obj.fun;
            end
            uMeshFun = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(f);
            obj.unfittedMeshFunction = uMeshFun;
        end

    end
end