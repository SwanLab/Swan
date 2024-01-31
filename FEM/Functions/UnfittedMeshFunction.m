classdef UnfittedMeshFunction < handle

    properties (Access = private)
        unfittedMesh
        levelSet
        cutCells
    end

    properties (Access = private)
        subMeshQuad
        fCutValues
    end

    properties (Access = public)
        innerMeshFunction
        innerCutMeshFunction
        boundaryCutMeshFunction
        unfittedBoundaryMeshFunction
        backgroundFunction
    end

    methods (Access = public)
        function obj = UnfittedMeshFunction(cParams)
            obj.init(cParams);
            obj.computeSubMeshQuadrature();
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
                        x1    = iMesh.coord(:,1);
                        y1    = iMesh.coord(:,2);
                        z1    = obj.innerMeshFunction.fValues;
                        a1    = trisurf(iMesh.connec,x1,y1,z1);
                    end
                    hold on
                    if ~isempty(obj.unfittedMesh.innerCutMesh)
                        iCMesh = obj.unfittedMesh.innerCutMesh.mesh;
                        x2     = iCMesh.coord(:,1);
                        y2     = iCMesh.coord(:,2);
                        z2     = obj.innerCutMeshFunction.fValues;
                        a2 = trisurf(iCMesh.connec,x2,y2,z2);
                    end
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
            obj.cutCells     = cParams.cutCells;
        end

        function computeSubMeshQuadrature(obj)
            mesh = obj.unfittedMesh.backgroundMesh;
            q    = Quadrature.set(mesh.type);
            q.computeQuadrature('CONSTANT');
            obj.subMeshQuad = q;
        end

        function computeBackgroundFunction(obj,f)
            fType = f.fType;
            switch fType
                case 'L2' % Provisional
                    f = f.project('P1');
            end
            obj.backgroundFunction = f;
        end

        function computeInnerMeshFunction(obj)
            if ~isempty(obj.unfittedMesh.innerMesh)
                fP1        = obj.backgroundFunction;
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
                cMeshGlobal              = obj.computeNonCutMesh();
                [fCutMesh,lsCutMesh]     = obj.computeNodalValuesFromNonCutMesh();
                innerValues              = fCutMesh(obj.isInterior(lsCutMesh));
                subMesh                  = obj.computeSubMesh(cMeshGlobal);
                subParams                = obj.computeNewNodalValuesFromSubMesh(fCutMesh,lsCutMesh,cMeshGlobal);
                subLevelSet              = [lsCutMesh;subParams.subMeshLevelSet];
                subfValues               = [fCutMesh;subParams.subMeshfValues];
                cutValues                = obj.computeCutValues(subMesh,subLevelSet,subfValues);
                s.mesh                   = obj.unfittedMesh.innerCutMesh.mesh;
                s.fValues                = [innerValues;subParams.cutMeshfValues;cutValues];
                fP1InnerCut              = P1Function(s);
                obj.innerCutMeshFunction = fP1InnerCut;
                obj.fCutValues           = cutValues;
            end
        end

        function computeBoundaryCutMeshFunction(obj)
            if ~isempty(obj.unfittedMesh.boundaryCutMesh)
                bCutMesh    = obj.unfittedMesh.boundaryCutMesh.mesh;
                ls          = obj.levelSet;
                f           = obj.backgroundFunction.fValues;
                innerValues = f(ls==0);
                cutValues   = obj.fCutValues;
                s.mesh      = bCutMesh;
                s.fValues   = [innerValues;cutValues];
                bCutFun     = P1Function(s);
                obj.boundaryCutMeshFunction = bCutFun;
            end
        end

        function computeUnfittedBoundaryMeshFunction(obj)
            uBoundMesh = obj.unfittedMesh.unfittedBoundaryMesh;
            fP1        = obj.backgroundFunction;
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

        function cMesh = computeNonCutMesh(obj)
            mesh     = obj.unfittedMesh.backgroundMesh;
            s.coord  = mesh.coord;
            s.connec = mesh.connec(obj.cutCells,:);
            s.kFace  = mesh.kFace;
            cMesh    = Mesh(s);
        end

        function [fCutMesh,lsCutMesh] = computeNodalValuesFromNonCutMesh(obj)
            fP1       = obj.backgroundFunction;
            mesh      = obj.unfittedMesh.backgroundMesh;
            nodes     = unique(mesh.connec(obj.cutCells,:));
            lsCutMesh = obj.levelSet(nodes);
            fCutMesh  = fP1.fValues(nodes);
        end

        function subMesh = computeSubMesh(obj,cMeshGlobal)
                switch obj.unfittedMesh.backgroundMesh.type
                    case 'QUAD'
                        ssub.mesh     = cMeshGlobal;
                        ssub.lastNode = obj.unfittedMesh.backgroundMesh.nnodes;
                        subMesher     = SubMesher(ssub);
                        subMesh       = subMesher.subMesh.computeCanonicalMesh();
                    case 'HEXAHEDRA'
                        Xiso     =  [-1 ,-1, -1;...
                            1, -1, -1;...
                            1, 1, -1;...
                            -1, 1, -1;...
                            -1, -1, 1;...
                            1, -1, 1;...
                            1, 1, 1;...
                            -1, 1, 1;];
                        connecIso    = delaunay(Xiso);
                        ss.coord     = Xiso;
                        ss.connec    = connecIso;
                        localMesh    = Mesh(ss);
                        nelem        = cMeshGlobal.nelem;
                        bCutConnec   = cMeshGlobal.connec;
                        connecIso    = localMesh.connec;
                        nElemIso     = size(connecIso,1);
                        nnodeSubMesh = size(connecIso,2);
                        subConnec    = bCutConnec(:,connecIso');
                        subConnec    = reshape(subConnec',[nnodeSubMesh,nelem*nElemIso])';
                        sss.coord    = cMeshGlobal.coord;
                        sss.connec   = subConnec;
                        subMesh      = Mesh(sss);
                        subMesh      = subMesh.computeCanonicalMesh();

                    otherwise
                        subMesh = cMeshGlobal.computeCanonicalMesh();
                end
        end

        function s = computeNewNodalValuesFromSubMesh(obj,fCutMesh,lsCutMesh,cMeshGlobal)
            switch obj.unfittedMesh.backgroundMesh.type
                case 'QUAD'
                    cMesh           = cMeshGlobal.computeCanonicalMesh();
                    sf.fValues      = fCutMesh;
                    sf.mesh         = cMesh;
                    fbMesh          = P1Function(sf);
                    sl.fValues      = lsCutMesh;
                    sl.mesh         = cMesh;
                    blsFun          = P1Function(sl);
                    fsubMesh        = obj.computeSubMeshValues(fbMesh);
                    lssubMesh       = obj.computeSubMeshValues(blsFun);
                    subMeshfValues  = fsubMesh.values;
                    subMeshLevelSet = lssubMesh.values;
                    cutMeshfValues  = subMeshfValues(obj.isInterior(subMeshLevelSet));
                case 'HEXAHEDRA'
                    subMeshfValues  = [];
                    subMeshLevelSet = [];
                    cutMeshfValues  = [];
                otherwise
                    subMeshfValues  = [];
                    subMeshLevelSet = [];
                    cutMeshfValues  = [];
            end
            s.subMeshfValues  = subMeshfValues;
            s.subMeshLevelSet = subMeshLevelSet;
            s.cutMeshfValues  = cutMeshfValues;
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

        function cutValues = computeCutValues(mesh,levelSet,fValues)
            m = mesh;
            m.computeEdges();
            e              = m.edges;
            s.nodesInEdges = e.nodesInEdges;
            s.levelSet     = levelSet;
            s.fValues      = fValues;
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