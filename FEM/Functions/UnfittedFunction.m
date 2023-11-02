classdef UnfittedFunction < handle

    properties (Access = public)
        ndimf
        unfittedMesh
        fun
    end

    methods (Access = public)

        function obj = UnfittedFunction(cParams)
            obj.init(cParams);
        end

        function fxV = evaluate(obj,xV)
            fxV = obj.evaluateInner(xV);
        end

        function fxV = evaluateInnerCut(obj,xV)
            inMesh = obj.unfittedMesh.innerCutMesh;
            obj.fun.updateMesh(inMesh.mesh);
            fxV    = obj.fun.evaluate(xV);
            %obj.fun.updateMesh(obj.unfittedMesh.backgroundMesh);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.unfittedMesh   = cParams.uMesh;
            obj.fun            = cParams.fun;
            obj.ndimf          = cParams.fun.ndimf;
        end

        function fxV = evaluateInner(obj,xV)
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