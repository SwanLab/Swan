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
                    file = 'PF_SENshear0_0025';
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
            end
            bc = phaseFieldBoundaryCreator(mesh,cParams);
        end

    end

end
