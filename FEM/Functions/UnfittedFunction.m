classdef UnfittedFunction < L2Function

    properties (Access = public)
        ndimf
        nDofs
        type
        unfittedMesh
    end

    properties (Access = private)
        fun
    end

    properties (Access = private)
        unfittedMeshFunction
    end

    methods (Access = public)

        function obj = UnfittedFunction(cParams)
            obj.init(cParams);
            obj.computeUnfittedMeshFunction();
        end

        function fxV = evaluate(obj,xV)
            fxV = obj.evaluateInnerElements(xV);
        end

        function fxV = evaluateCutElements(obj,q)
            uMeshFun = obj.unfittedMeshFunction;
            fNew     = uMeshFun.innerCutMeshFunction;
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
            uMeshFun = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(obj.fun);
            obj.unfittedMeshFunction = uMeshFun;
            obj.nDofs                = uMeshFun.backgroundFunction.nDofs;
            obj.type                 = uMeshFun.backgroundFunction.type;
        end

        function fxV = evaluateInnerElements(obj,xV)
            fxV     = obj.fun.evaluate(xV);
            gMesh   = obj.unfittedMesh.backgroundMesh;
            inMesh  = obj.unfittedMesh.innerMesh;
            gCoor   = gMesh.computeXgauss(xV);
            inCoor  = inMesh.mesh.computeXgauss(xV);
            gCoor1  = squeeze(gCoor(:,1,:))';
            inCoor1 = squeeze(inCoor(:,1,:))';
            isVoid  = not(ismember(gCoor1,inCoor1,'rows'));
            fxV(:,:,isVoid) = 0;
        end

    end
end