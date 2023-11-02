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
            mesh      = obj.unfittedMesh.backgroundMesh;
            bCMesh    = obj.unfittedMesh.boundaryCutMesh;
            connec    = mesh.connec;
            inCConnec = connec(bCMesh.cellContainingSubcell,:);
            s.connec  = inCConnec;
            s.coord   = bCMesh.mesh.coord;
            meshNew   = Mesh(s);
            obj.fun.updateMesh(meshNew);
            fxV       = obj.fun.evaluate(xV);
            obj.fun.updateMesh(mesh);
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