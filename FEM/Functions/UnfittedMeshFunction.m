classdef UnfittedMeshFunction < handle

    properties (Access = private)
        unfittedMesh
        levelSet
    end

    properties (Access = private)
        subMeshQuad
        funP1
    end

    properties (Access = public)
        innerMeshFunction
        innerCutMeshFunction
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
        end

        function plot(obj)
            switch obj.unfittedMesh.backgroundMesh.type
                case {'TRIANGLE','QUAD'}
                    iMesh  = obj.unfittedMesh.innerMesh.mesh;
                    iCMesh = obj.unfittedMesh.innerCutMesh.mesh;
                    x1     = iMesh.coord(:,1);
                    y1     = iMesh.coord(:,2);
                    x2     = iCMesh.coord(:,1);
                    y2     = iCMesh.coord(:,2);
                    z1     = obj.innerMeshFunction.fValues;
                    z2     = obj.innerCutMeshFunction.fValues;
                    figure()
                    %for idim = 1:obj.ndimf
                        %subplot(1,obj.ndimf,idim);
                        a1 = trisurf(iMesh.connec,x1,y1,z1);
                        hold on
                        a2 = trisurf(iCMesh.connec,x2,y2,z2);
                        hold off
                        view(0,90)
                        shading interp
                        a1.EdgeColor = [0 0 0];
                        a2.EdgeColor = [0 0 0];
                    %end

                otherwise

            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.unfittedMesh = cParams.uMesh;
            obj.levelSet     = cParams.levelSet;
        end

        function computeSubMeshQuadrature(obj)
            mesh = obj.unfittedMesh.backgroundMesh;
            q    = Quadrature.set(mesh.type);
            q.computeQuadrature('CONSTANT');
            obj.subMeshQuad = q;
        end

        function computeInnerMeshFunction(obj)
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

        function computeInnerCutMeshFunction(obj)
            switch obj.unfittedMesh.backgroundMesh.type
                case 'QUAD'
                    obj.computeCutMeshFunctionQuadrilateral();
                case 'HEXAHEDRA'
                    obj.computeCutMeshFunctionHexahedra();
                otherwise
                    obj.computeCutMeshFunctionStandard();
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
        end

        % hexahedra

        function computeCutMeshFunctionStandard(obj)
            iCMesh           = obj.unfittedMesh.innerCutMesh;
            subCells         = unique(iCMesh.cellContainingSubcell);
            innerValues      = obj.computeInnerValuesFromCutMesh(subCells);
            c.mesh           = obj.unfittedMesh.backgroundMesh;
            c.levelSet       = obj.levelSet;
            c.fValues        = obj.funP1.fValues;
            cutValues        = obj.computeCutValues(c);
            % ME QUEDÉ AQUÍ
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