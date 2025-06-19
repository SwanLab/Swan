classdef BenchmarkManager < handle

    methods (Access = public, Static)

        function [mesh, bc] = create(cParams)
            switch cParams.type.mesh
                case '1Elem'
                    mesh = QuadMesh(1,1,1,1);
                case 'nElem'
                    N = cParams.N;
                    mesh = QuadMesh(1,1,N,N);
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
                    file = 'Hole0_0025';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
            end
            bc = PhaseFieldBoundaryCreator(mesh,cParams);
        end

    end

end
