classdef CutMeshFunctionComputer < handle

    properties (Access = private)
        backgroundMesh
        cutMesh
    end

    properties (Access = private)
        coordGlobal
        edgesBackgr
        edgesCutGlobal
        parallelEdges
    end

    methods (Access = public)
        function obj = CutMeshFunctionComputer(cParams)
            obj.init(cParams);
            obj.computeGlobalCoordinatesMatrix();
            obj.computeBackgroundMeshEdges();
            obj.computeCutMeshEdgesGlobal();
            obj.computeParallelEdgesMatrix();
        end

        function newFeFun = compute(obj,backgroundFeFun)
            fB               = backgroundFeFun;
            oldfValuesCMesh  = obj.computeRemainingOldFValuesAtCutMesh(fB);
            addfValues       = obj.computeNewFValues(fB);
            newfValues       = [oldfValuesCMesh;addfValues];
            s.feFunType      = class(fB);
            s.mesh           = obj.cutMesh;
            s.ndimf          = fB.ndimf;
            newFeFun         = FeFunction.createEmpty(s);
            newFeFun.fValues = newfValues;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.cutMesh        = cParams.cutMesh;
        end

        function computeGlobalCoordinatesMatrix(obj)
            coord           = obj.backgroundMesh.coord;
            icCoord         = obj.cutMesh.coord;
            ind             = not(ismember(icCoord,coord,'rows'));
            obj.coordGlobal = [coord;icCoord(ind,:)];
        end

        function computeBackgroundMeshEdges(obj)
            mesh            = obj.backgroundMesh;
            eBackgr         = obj.computeEdgesSorted(mesh);
            obj.edgesBackgr = eBackgr;
        end

        function computeCutMeshEdgesGlobal(obj)
            mesh         = obj.backgroundMesh;
            cMesh        = obj.cutMesh;
            e            = obj.computeEdgesSorted(cMesh);
            cBackgr      = mesh.coord;
            cCMesh       = cMesh.coord;
            ind          = find(not(ismember(cCMesh,cBackgr,'rows')));
            ind2         = find(not(ismember(cBackgr,cCMesh,'rows')));
            n            = numel(ind2);
            e            = obj.keepCutMeshEdgesCoincidentToBackgrEdges(e);
            e(e>=ind(1)) = e(e>=ind(1))+n;
            obj.edgesCutGlobal = e;
        end

        function e = keepCutMeshEdgesCoincidentToBackgrEdges(obj,e)
            cBackgr = obj.backgroundMesh.coord;
            cCMesh  = obj.cutMesh.coord;
            ind     = find(not(ismember(cCMesh,cBackgr,'rows')));
            ind     = reshape(ind,[1,1,length(ind)]);
            b       = sum(e==ind,3);
            e       = e(xor(b(:,1),b(:,2)),:);
        end

        function computeParallelEdgesMatrix(obj)
            eB                = obj.edgesBackgr;
            e                 = obj.edgesCutGlobal;
            vB                = obj.computeEdgesVectorMatrix(eB);
            v                 = obj.computeEdgesVectorMatrix(e);
            parallelBgrNodes  = (abs(v*vB'-1)<=1e-6)*eB;
            eFinal            = [e,parallelBgrNodes(:,2)];
            notZeros          = not(sum(ismember(eFinal,0),2));
            eFinal            = eFinal(notZeros,:);
            obj.parallelEdges = eFinal;
        end

        function v = computeEdgesVectorMatrix(obj,e)
            coordGlob = obj.coordGlobal;
            e1        = e(:,1);
            e2        = e(:,2);
            v         = coordGlob(e2,:)-coordGlob(e1,:);
            v         = v./sqrt(sum(v.^2,2));
        end

        function f = computeRemainingOldFValuesAtCutMesh(obj,fB)
            oldfValues = fB.fValues;
            cBackgr    = obj.backgroundMesh.coord;
            cCMesh     = obj.cutMesh.coord;
            ind        = ismember(cBackgr,cCMesh,'rows');
            f          = oldfValues(ind,:);
        end

        function f = computeNewFValues(obj,fB)
            coordGlob = obj.coordGlobal;
            eFinal    = obj.parallelEdges;
            nodes1    = eFinal(:,1);
            nodes2    = eFinal(:,2);
            nodes3    = eFinal(:,3);
            xi        = coordGlob(nodes1,:);
            x         = coordGlob(nodes2,:);
            xf        = coordGlob(nodes3,:);
            l1        = sqrt(sum((x-xi).^2,2));
            l2        = sqrt(sum((xf-x).^2,2));
            f1        = fB.fValues(nodes1,:);
            f2        = fB.fValues(nodes3,:);
            fNew      = f2+(f1-f2).*(l2./(l1+l2));
            [~,v]     = unique(nodes2);
            f         = fNew(v,:);
        end
    end

    methods (Access = private, Static)
        function e = computeEdgesSorted(m)
            m.computeEdges();
            e                  = m.edges.nodesInEdges;
            V                  = e;
            e(V(:,1)>V(:,2),1) = V(V(:,1)>V(:,2),2);
            e(V(:,1)>V(:,2),2) = V(V(:,1)>V(:,2),1);
        end
    end
end