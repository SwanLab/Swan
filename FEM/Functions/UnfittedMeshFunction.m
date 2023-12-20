classdef UnfittedMeshFunction < handle

    properties (Access = private)
        unfittedMesh
        levelSet
        cutCells
    end

    properties (Access = private)
        subMeshQuad
        funP1
        fCutValues
    end

    properties (Access = public)
        innerMeshFunction
        innerCutMeshFunction
        boundaryCutMeshFunction
        unfittedBoundaryMeshFunction
    end

    methods (Access = public)
        function obj = UnfittedMeshFunction(cParams)
            obj.init(cParams);
            obj.computeSubMeshQuadrature();
        end

        function compute(obj,f)
            obj.funP1 = f.project('P1');
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
                        x1    = iMesh.coord(:,1);
                        y1    = iMesh.coord(:,2);
                        z1    = obj.innerMeshFunction.fValues;
                        a1    = trisurf(iMesh.connec,x1,y1,z1);
                    end
                    hold on
                    iCMesh = obj.unfittedMesh.innerCutMesh.mesh;
                    x2     = iCMesh.coord(:,1);
                    y2     = iCMesh.coord(:,2);
                    z2     = obj.innerCutMeshFunction.fValues;
                    a2 = trisurf(iCMesh.connec,x2,y2,z2);
                    hold off
                    view(0,90)
                    shading interp
                    a1.EdgeColor = [0 0 0];
                    a2.EdgeColor = [0 0 0];

                otherwise

            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.unfittedMesh = cParams.uMesh;
            obj.levelSet     = cParams.levelSet;
            if isfield(cParams,'cutCells')
                obj.cutCells = cParams.cutCells;
            end
        end

        function computeSubMeshQuadrature(obj)
            mesh = obj.unfittedMesh.backgroundMesh;
            q    = Quadrature.set(mesh.type);
            q.computeQuadrature('CONSTANT');
            obj.subMeshQuad = q;
        end

        function computeInnerMeshFunction(obj)
            if ~isempty(obj.unfittedMesh.innerMesh)
                fP1        = obj.funP1;
                iMesh      = obj.unfittedMesh.innerMesh;
                connecLoc  = iMesh.mesh.connec;
                connecGlob = iMesh.globalConnec;
                glob2loc(connecLoc(:)) = connecGlob(:);
                s.fValues = fP1.fValues(glob2loc);
                s.mesh    = iMesh.mesh;
                fP1Inner  = P1Function(s);
                obj.innerMeshFunction = fP1Inner;
            end
        end

        function computeInnerCutMeshFunction(obj)
            if ~isempty(obj.unfittedMesh.innerCutMesh)
                switch obj.unfittedMesh.backgroundMesh.type
                    case 'QUAD'
                        obj.computeCutMeshFunctionQuadrilateral();
                    case 'HEXAHEDRA'
                        obj.computeCutMeshFunctionStandard();
                    otherwise
                        obj.computeCutMeshFunctionStandard();
                end
            end
        end

        function computeBoundaryCutMeshFunction(obj)
            if ~isempty(obj.unfittedMesh.boundaryCutMesh)
                bCutMesh    = obj.unfittedMesh.boundaryCutMesh.mesh;
                ls          = obj.levelSet;
                f           = obj.funP1.fValues;
                innerValues = f(ls==0);
                cutValues   = obj.fCutValues;
                s.mesh      = bCutMesh;
                s.fValues   = [innerValues;cutValues];
                bCutFun     = P1Function(s);
                obj.boundaryCutMeshFunction = bCutFun;
            end
        end

        function computeUnfittedBoundaryMeshFunction(obj)
            uBoundMesh   = obj.unfittedMesh.unfittedBoundaryMesh;
            fP1          = obj.funP1;
            if ~isempty(uBoundMesh.meshes)
                activeMeshes = uBoundMesh.getActiveMesh();
                for i = 1:length(activeMeshes)
                    uMeshi     = activeMeshes{i};
                    connecLoc  = uMeshi.backgroundMesh.connec;
                    connecGlob = uBoundMesh.getGlobalConnec{i};
                    glob2loc(connecLoc(:)) = connecGlob(:);
                    s.fValues  = fP1.fValues(glob2loc);
                    s.mesh     = uMeshi.backgroundMesh;
                    fbackMeshi = P1Function(s);
                    glob2loc   = [];
                    fi         = uMeshi.obtainFunctionAtUnfittedMesh(fbackMeshi);
                    obj.unfittedBoundaryMeshFunction.activeFuns{i} = fi;
                end
            end
        end

        function computeCutMeshFunctionQuadrilateral(obj)
            mesh             = obj.unfittedMesh.backgroundMesh;
            sls.fValues      = obj.levelSet;
            sls.mesh         = mesh;
            fls              = P1Function(sls);
            iCMesh           = obj.unfittedMesh.innerCutMesh;
            subCells         = unique(iCMesh.cellContainingSubcell);
            innerValues      = obj.computeInnerValuesFromCutMesh(subCells);
            lsSubMesh        = obj.computeSubMeshValues(fls);
            subMesh          = obj.computeSubMeshValues(obj.funP1);
            lsSubCutMesh     = lsSubMesh.values(subCells);
            subCutMeshValues = subMesh.values(subCells);
            subCutMeshValues = subCutMeshValues(obj.isInterior(lsSubCutMesh));
            ssub.mesh        = mesh;
            ssub.lastNode    = mesh.nnodes;
            subMesher        = SubMesher(ssub);
            c.mesh           = subMesher.subMesh;
            c.levelSet       = lsSubMesh.allValues;
            c.fValues        = subMesh.allValues;
            cutValues        = obj.computeCutValues(c);
            fValues          = [innerValues;subCutMeshValues;cutValues];
            ss.mesh          = iCMesh.mesh;
            ss.fValues       = fValues;
            fP1InnerCut      = P1Function(ss);
            obj.innerCutMeshFunction = fP1InnerCut;
            obj.fCutValues           = cutValues;
        end

        function computeCutMeshFunctionHexahedra(obj)
            cutPointsCalculator   = CutPointsCalculator();
            s.backgroundCutCells  = obj.cutCells;
            s.backgroundMesh      = obj.unfittedMesh.backgroundMesh;
            s.levelSet_background = obj.levelSet;
            cutPointsCalculator.init(s);
            cutPointsCalculator.computeCutPoints();
            connec    = obj.unfittedMesh.backgroundMesh.connec;
            fValues   = [];
            cutValues = [];
            coorGlob  = [];
            cutCoord  = [];
            for i = 1:length(obj.cutCells)
                nodes     = connec(obj.cutCells(i),:)';
                isActive  = obj.isInterior(obj.levelSet(nodes));
                dofs      = nodes(isActive);
                xV        = cutPointsCalculator.getThisCellCutPoints(i).ISO';
                fxV       = obj.funP1.evaluate(xV);
                fxV       = fxV(:,:,obj.cutCells(i))';
                xxV       = cutPointsCalculator.getThisCellCutPoints(i).GLOBAL;
                fValues   = [fValues;obj.funP1.fValues(dofs,:);fxV];
                cutValues = [cutValues;fxV];
                coorGlob  = [coorGlob;obj.unfittedMesh.backgroundMesh.coord(dofs,:);xxV];
                cutCoord  = [cutCoord;xxV];
            end
            [~,v]       = unique(coorGlob,'stable','rows');
            fValues     = fValues(v);
            [~,v]       = unique(cutCoord,'stable','rows');
            cutValues   = cutValues(v);
            ss.mesh     = obj.unfittedMesh.innerCutMesh.mesh;
            ss.fValues  = fValues;
            fP1InnerCut = P1Function(ss);
            obj.innerCutMeshFunction = fP1InnerCut;
            obj.fCutValues           = cutValues;
        end

        function computeCutMeshFunctionStandard(obj)
            iCMesh           = obj.unfittedMesh.innerCutMesh;
            subCells         = unique(iCMesh.cellContainingSubcell);
            innerValues      = obj.computeInnerValuesFromCutMesh(subCells);
            c.mesh           = obj.unfittedMesh.backgroundMesh;
            c.levelSet       = obj.levelSet;
            c.fValues        = obj.funP1.fValues;
            cutValues        = obj.computeCutValues(c);
            fValues          = [innerValues;cutValues];
            ss.mesh          = iCMesh.mesh;
            ss.fValues       = fValues;
            fP1InnerCut      = P1Function(ss);
            obj.innerCutMeshFunction = fP1InnerCut;
            obj.fCutValues           = cutValues;
        end

        function innerValues = computeInnerValuesFromCutMesh(obj,subCells)
            fP1         = obj.funP1;
            mesh        = obj.unfittedMesh.backgroundMesh;
            nodes       = unique(mesh.connec(subCells,:));
            ls          = obj.levelSet(nodes);
            innerNodes  = nodes(obj.isInterior(ls));
            innerValues = fP1.fValues(innerNodes);
        end

        function s = computeSubMeshValues(obj,fP1)
            q           = obj.subMeshQuad;
            xV          = q.posgp;
            s.values    = squeeze(fP1.evaluate(xV));
            s.allValues = [fP1.fValues;s.values];
        end

    end

    methods (Static, Access = private)

        function isNodeInterior = isInterior(ls)
            isNodeInterior = logical(1-heaviside(ls));
        end

        function cutValues = computeCutValues(cParams)
            m = cParams.mesh;
            m.computeEdges();
            e              = m.edges;
            s.nodesInEdges = e.nodesInEdges;
            s.levelSet     = cParams.levelSet;
            s.fValues      = cParams.fValues;
            ce             = CutEdgesComputer(s);
            ce.compute();
            s.xCutEdgePoint = ce.xCutEdgePoint;
            s.isEdgeCut     = ce.isEdgeCut;
            cf              = CutFunctionValuesComputer(s);
            cf.compute();
            cutValues = cf.cutValues;
        end
    end
end