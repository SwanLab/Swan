classdef BenchmarkManager < handle

    methods (Access = public, Static)

        function [mesh, bc] = create(cParams)
            switch cParams.type.mesh
                case '1Elem'
                    mesh = QuadMesh(1,1,1,1);
                case 'nElem'
                    N = cParams.N;
                    mesh = QuadMesh(1,10,N,N);
                case 'SENshear'
                    file = 'PF_SENshear0_01';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'SENtraction'
                    file = 'PF_SENtraction0_0025';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'SENmixed'
                    file = 'PF_SENmixed0_0025';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'Lshape'
                    file = 'PF_Lmesh';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'FiberMatrix'
                    bgMesh = TriangleMesh(1,1,20,20);
                    sLS.type        = 'CircleInclusion';
                    sLS.xCoorCenter = 0.5;
                    sLS.yCoorCenter = 0.5;
                    sLS.radius      = 0.2;
                    g               = GeometricalFunction(sLS);
                    lsFun           = g.computeLevelSetFunction(bgMesh);
                    levelSet        = lsFun.fValues;
                    sUm.backgroundMesh = bgMesh;
                    sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
                    uMesh = UnfittedMesh(sUm);
                    uMesh.compute(levelSet);
                    mesh = uMesh.createInnerMeshGoodConditioning();
            end
            bc = phaseFieldBoundaryCreator(mesh,cParams);
        end

    end

end
