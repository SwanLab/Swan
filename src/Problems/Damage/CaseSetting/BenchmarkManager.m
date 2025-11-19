classdef BenchmarkManager < handle

    methods (Access = public, Static)

        function [mesh, bc] = create(cParams)
            switch cParams.mesh.type
                case '1Elem'
                    mesh = QuadMesh(1,1,1,1);
                case 'Rectangle'
                    l = cParams.mesh.length;
                    w = cParams.mesh.width;
                    N = cParams.mesh.lN;
                    M = cParams.mesh.wN;
                    mesh = QuadMesh(l,w,N,M);
                case 'SENtest'
                    file = 'SENtest0_05';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'Sample'
                    file = 'Sample';
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
            end
            bc.u   = BoundaryConditionsCreator(mesh,cParams.bc.u);
            bc.phi = BoundaryConditionsCreator(mesh,cParams.bc.phi);
        end

    end

end
