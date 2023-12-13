classdef UnfittedMeshFunction < handle

    properties (Access = private)
        unfittedMesh
        levelSet
    end

    properties (Access = public)
        innerMeshFunction
        innerCutMeshFunction
    end

    methods (Access = public)
        function obj = UnfittedMeshFunction(cParams)
            obj.init(cParams);
        end

        function compute(obj,f)
            fP1 = f.project('P1');
            obj.computeInnerMeshFunction(fP1);
            obj.computeInnerCutMeshFunction(fP1);
        end

        function plot(obj)
            switch obj.unfittedMesh.backgroundMesh.type
                case {'TRIANGLE','QUAD'}
                    x1 = obj.unfittedMesh.innerMesh.mesh.coord(:,1);
                    y1 = obj.unfittedMesh.innerMesh.mesh.coord(:,2);
                    x2 = obj.unfittedMesh.innerCutMesh.mesh.coord(:,1);
                    y2 = obj.unfittedMesh.innerCutMesh.mesh.coord(:,2);
                    x  = [x1;x2];
                    y  = [y1;y2];
                    figure()
                    %for idim = 1:obj.ndimf
                        subplot(1,obj.ndimf,idim);
                        z = obj.fValues(:,idim);
                        a = trisurf(obj.mesh.connec,x,y,z);
                        view(0,90)
                        %             colorbar
                        shading interp
                        a.EdgeColor = [0 0 0];
                        title(['dim = ', num2str(idim)]);
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

        function computeInnerMeshFunction(obj,fP1)
            iMesh      = obj.unfittedMesh.innerMesh;
            connecLoc  = iMesh.mesh.connec;
            connecGlob = iMesh.globalConnec;
            glob2loc(connecLoc(:)) = connecGlob(:);
            s.fValues = fP1.fValues(glob2loc);
            s.mesh    = iMesh.mesh;
            fP1Inner  = P1Function(s);
            obj.innerMeshFunction = fP1Inner;
        end

        function computeInnerCutMeshFunction(obj,fP1)
            switch obj.unfittedMesh.backgroundMesh.type
                case 'QUAD'
                    obj.computeCutMeshFunctionQuadrilateral(fP1);
                case 'HEXAHEDRA'
                    obj.computeCutMeshFunctionHexahedra();
                otherwise
                    obj.computeCutMeshFunctionProvisional();
            end
        end

        function computeCutMeshFunctionQuadrilateral(obj,fP1)
            iCMesh         = obj.unfittedMesh.innerCutMesh;
            mesh           = obj.unfittedMesh.backgroundMesh;
            subCells       = unique(iCMesh.cellContainingSubcell);
            nodes          = unique(mesh.connec(subCells,:));
            ls             = obj.levelSet(nodes);
            isNodeInterior = logical(1-heaviside(ls));
            innerValues    = fP1.fValues(nodes(isNodeInterior));


            q = Quadrature.set(mesh.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            sls.fValues = obj.levelSet;
            sls.mesh    = mesh;
            fls         = P1Function(sls);
            lsSubMesh   = squeeze(fls.evaluate(xV));
            lsSubCutMesh   = lsSubMesh(subCells);
            subMeshValues = squeeze(fP1.evaluate(xV));
            subCutMeshValues = subMeshValues(subCells);
            isNodeInterior = logical(1-heaviside(lsSubCutMesh));
            subCutMeshValues = subCutMeshValues(isNodeInterior);


            ssub.mesh        = mesh;
            ssub.lastNode    = mesh.nnodes;
            subMesher = SubMesher(ssub);
            subMesher.subMesh.computeEdges();
            e = subMesher.subMesh.edges;
            s.nodesInEdges = e.nodesInEdges;
            s.levelSet     = [obj.levelSet;lsSubMesh];
            s.fValues       = [fP1.fValues;subMeshValues];


            ce             = CutEdgesComputer(s);
            ce.compute();
            s.xCutEdgePoint = ce.xCutEdgePoint;
            s.isEdgeCut     = ce.isEdgeCut;
            cf              = CutFunctionValuesComputer(s);
            cf.compute();
            fValues = cf.cutValues;
            fValues              = [innerValues;subCutMeshValues;fValues];


            ss.mesh = iCMesh.mesh;
            ss.fValues = fValues;
            fP1InnerCut = P1Function(ss);
            obj.innerCutMeshFunction = fP1InnerCut;
        end
    end
end