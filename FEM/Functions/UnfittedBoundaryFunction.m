classdef UnfittedBoundaryFunction < handle

    properties (Access = public)
        ndimf
        unfittedMesh
        fun
    end

    methods (Access = public)

        function obj = UnfittedBoundaryFunction(cParams)
            obj.init(cParams);
        end

        function fxV = evaluate(obj,xV)
            % Will be analogous to inner in UnfittedFunction and then the whole boundaryUnfitted will be necessary
        end

        function fxV = evaluateCutElements(obj,xV)
            % FeFunction:
            mesh          = obj.unfittedMesh.backgroundMesh;
            coordOriginal = mesh.coord;
            n0            = size(coordOriginal,1)+1;
            bCMesh        = obj.unfittedMesh.boundaryCutMesh;
            coordComplete = bCMesh.mesh.coord;
            oldfValues = obj.fun.fValues;
            x          = coordOriginal(:,1);
            y          = coordOriginal(:,2);
            F          = scatteredInterpolant(x,y,oldfValues);
            newfValues = F(coordComplete(n0:end,:));
            s.fValues = [oldfValues;newfValues];
            s.mesh    = bCMesh.mesh;
            newFun    = eval([class(obj.fun),'(s)']);
            fxV       = newFun.evaluate(xV);

            % L2:
%             obj.fun.updateMesh(bCMesh.mesh);
%             fxV       = obj.fun.evaluate(xV);
%             obj.fun.updateMesh(mesh);
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