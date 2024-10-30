classdef benchmarkManager < handle

    methods (Access = public, Static)

        function [mesh, bc] = create(cParams)
            switch cParams.type.mesh
                case '1Elem'
                    mesh = QuadMesh(1,1,1,1);
                case 'nElem'
                    N = cParams.N;
                    mesh = QuadMesh(1,10,N,N);
                case 'SEN'
                    file = 'SENmeshDisplaced';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'Lshape'
                    file = 'PF_Lmesh';
                    a.fileName = file;
                    s = FemDataContainer(a);
                    mesh = s.mesh;
                case 'FiberMatrix'
                    % %%%%%% REVIEW LEVEL SETS %%%%%%%%%%%%%%%
                    % % Generate coordinates
                    % x1 = linspace(0,1,20);
                    % x2 = linspace(1,2,20);
                    % % Create the grid
                    % [xv,yv] = meshgrid(x1,x2);
                    % % Triangulate the mesh to obtain coordinates and connectivities
                    % [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
                    % sBg.coord  = V(:,1:2);
                    % sBg.connec = F;
                    % bgMesh = Mesh.create(sBg);
                    % bdMesh  = bgMesh.createBoundaryMesh();
                    %
                    % % Level set creation
                    % sLS.type       = 'circleInclusion';
                    % sLS.mesh       = bgMesh;
                    % sLS.ndim       = 2;
                    % sLS.fracRadius = 0.4;
                    % sLS.coord      = bgMesh.coord;
                    % ls = LevelSetCreator.create(sLS);
                    % levelSet = ls.getValue();
                    %
                    % sUm.backgroundMesh = bgMesh;
                    % sUm.boundaryMesh   = bdMesh;
                    % uMesh = UnfittedMesh(sUm);
                    % uMesh.compute(levelSet);
                    %
                    % mesh = uMesh.createInnerMeshGoodConditioning();
            end
            bc = phaseFieldBoundaryCreator(mesh,cParams);
        end

    end

end
