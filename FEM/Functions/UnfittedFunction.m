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

        function fxV = evaluateCutElements(obj,xV)
            mesh      = obj.unfittedMesh.backgroundMesh;
            inCMesh   = obj.unfittedMesh.innerCutMesh;
            connec    = mesh.connec;
            if isempty(inCMesh)
                fxV = zeros(1,100,1);
            else
                inCConnec = connec(inCMesh.cellContainingSubcell,:);
                s.connec  = inCConnec;
                s.coord   = inCMesh.mesh.coord;
                if size(s.coord,2)==2 && size(s.connec,2)==2
                    s.kFace = -1;
                end
                if size(s.coord,2)==3 && size(s.connec,2)==3
                    s.kFace = -1;
                end
                meshNew   = Mesh(s);
                obj.fun.updateMesh(meshNew);
                fxV       = obj.fun.evaluate(xV);
                obj.fun.updateMesh(mesh);
            end
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
            if isempty(inMesh)
                fxV(:,:,:) = 0;
            else
                gCoor   = gMesh.computeXgauss(xV);
                inCoor  = inMesh.mesh.computeXgauss(xV);
                gCoor1  = squeeze(gCoor(:,1,:))';
                inCoor1 = squeeze(inCoor(:,1,:))';
                isVoid  = not(ismember(gCoor1,inCoor1,'rows'));
                fxV(:,:,isVoid) = 0;
            end
        end

    end
end