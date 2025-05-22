classdef BenchmarkManager < handle

    methods (Access = public, Static)

        function [mesh, bc] = create(cParams)
            switch cParams.type.mesh
                case '1Elem'
                    mesh = QuadMesh(1,1,1,1);
                case 'nElem'
                    N = cParams.N;
                    mesh = QuadMesh(1,10,N,N);
                case 'SENtest'
                    file = 'SENtest0_05';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'SENshear'
                    file = 'SENshear0_0025';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'SENtraction'
                    file = 'SENtraction0_0025';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'SENmixed'
                    file = 'SENmixed0_0025';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case {'FiberMatrix','Hole'}
                    N = 100;
                    mesh = QuadMesh(1,1,N,N);

                    gPar.type = 'Circle';
                    gPar.xCoorCenter = 0.5;
                    gPar.yCoorCenter = 0.5;
                    gPar.radius = 0.25;
                    g                  = GeometricalFunction(gPar);
                    ls = g.computeLevelSetFunction(mesh);
                    sUm.backgroundMesh = mesh;
                    sUm.boundaryMesh   = mesh.createBoundaryMesh;
                    uMesh              = UnfittedMesh(sUm);
                    uMesh.compute(-ls.fValues);
                    mesh = uMesh.createInnerMeshGoodConditioning();
                    mesh.computeCanonicalMesh()
            end
            bc = PhaseFieldBoundaryCreator(mesh,cParams);
        end

    end

end
