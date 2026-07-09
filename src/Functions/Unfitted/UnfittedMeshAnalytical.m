classdef UnfittedMeshAnalytical < handle

    properties (Access = private)
        unfittedMesh
    end

    properties (Access = public)
        innerMeshFunction
        innerCutMeshFunction
        boundaryCutMeshFunction
        unfittedBoundaryMeshFunction
        backgroundFunction
    end

    methods (Access = public)
        function obj = UnfittedMeshAnalytical(cParams)
            obj.init(cParams);
        end

        function compute(obj,f)
            obj.computeBackgroundFunction(f);
            obj.computeInnerMeshFunction();
            obj.computeInnerCutMeshFunction();
            obj.computeBoundaryCutMeshFunction();
            obj.computeUnfittedBoundaryMeshFunction();
        end

        function plot(obj)
            switch obj.unfittedMesh.backgroundMesh.type
                case {'TRIANGLE','QUAD'}
                    figure()
                    if ~isempty(obj.unfittedMesh.innerMesh)
                        iMesh = obj.unfittedMesh.innerMesh.mesh;
                        iMF   = obj.innerMeshFunction.project('P1');
                        x1    = iMesh.coord(:,1);
                        y1    = iMesh.coord(:,2);
                        z1    = iMF.fValues;
                        a1    = trisurf(iMesh.connec,x1,y1,z1);
                    end
                    hold on
                    if ~isempty(obj.unfittedMesh.innerCutMesh)
                        iCMesh = obj.unfittedMesh.innerCutMesh.mesh;
                        iCMF   = obj.innerCutMeshFunction.project('P1');
                        x2     = iCMesh.coord(:,1);
                        y2     = iCMesh.coord(:,2);
                        z2     = iCMF.fValues;
                        a2 = trisurf(iCMesh.connec,x2,y2,z2);
                    end
                    hold off
                    view(0,90)
                    shading interp
                    a1.EdgeColor = [0 0 0];
                    a2.EdgeColor = [0 0 0];
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.unfittedMesh = cParams.uMesh;
        end

        function computeBackgroundFunction(obj,f)
            obj.backgroundFunction = f;
        end

        function computeInnerMeshFunction(obj)
            iM = obj.unfittedMesh.innerMesh;
            if ~isempty(iM)
                f      = obj.backgroundFunction;
                fInner = f.createNew(iM.mesh);
                obj.innerMeshFunction = fInner;
            end
        end

        function computeInnerCutMeshFunction(obj)
            iCM = obj.unfittedMesh.innerCutMesh;
            if ~isempty(iCM)
                f         = obj.backgroundFunction;
                fInnerCut = f.createNew(iCM.mesh);
                obj.innerCutMeshFunction = fInnerCut;
            end
        end

        function computeBoundaryCutMeshFunction(obj)
            bCM = obj.unfittedMesh.boundaryCutMesh;
            if ~isempty(bCM)
                f       = obj.backgroundFunction;
                bCutFun = f.createNew(bCM.mesh);
                obj.boundaryCutMeshFunction = bCutFun;
            end
        end

        function computeUnfittedBoundaryMeshFunction(obj)
            uF = [];
            f  = obj.backgroundFunction;
            uBoundMesh = obj.unfittedMesh.unfittedBoundaryMesh;
            if ~isempty(uBoundMesh.meshes)
                aMeshes = uBoundMesh.getActiveMesh();
                for i = 1:length(aMeshes)
                    uMeshi = aMeshes{i};
                    meshi  = uMeshi.backgroundMesh;
                    fi     = f.createNew(meshi);
                    uF.activeFuns{i} = uMeshi.obtainFunctionAtUnfittedMesh(fi);
                end
            end
            obj.unfittedBoundaryMeshFunction = uF;
        end
    end
end