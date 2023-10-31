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
            gMesh     = obj.unfittedMesh.backgroundMesh;
            inMesh    = obj.unfittedMesh.innerMesh;
            %inCutMesh = obj.unfittedMesh.innerCutMesh;

            gCoor   = gMesh.computeXgauss(xV);
            inCoor  = inMesh.mesh.computeXgauss(xV);
            gCoor1  = squeeze(gCoor(:,1,:))';
            inCoor1 = squeeze(inCoor(:,1,:))';
            isVoid = not(ismember(gCoor1,inCoor1,'rows'));

            fxV             = obj.fun.evaluate(xV);
            fxV(:,:,isVoid) = 0;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.unfittedMesh   = cParams.uMesh;
            obj.fun            = cParams.fun;
            obj.ndimf          = cParams.fun.ndimf;
        end

    end
end