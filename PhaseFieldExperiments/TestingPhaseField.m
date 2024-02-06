classdef TestingPhaseField < handle

    properties (Access = public)
        E = 210;
        nu = 0.3;
        Gc = 5e-3;
        l0 = 0.1;
        pExp = 2;
        % bcVal = [linspace(0,-2e-1,100), ...
        %         linspace(-2e-1,2e-1,400), ...
        %         linspace(2e-1,4e-1,100), ...
        %         linspace(4e-1,-4e-1,100), ...
        %         linspace(-4e-1,1,400)];
        bcVal = linspace(0,1,500);
        % bcVal = [0.001];
    end

    properties (Access = private)
    end

    properties (Access = private)
        mesh
        initialPhaseField
        materialPhaseField
        dissipationPhaseField
        Constant
    end

    methods (Access = public)

        function obj = TestingPhaseField()
            obj.init()
            obj.createMesh();
            obj.createInitialPhaseField();
            obj.createMaterialPhaseField();
            obj.createDissipationInterpolation();
            obj.solveProblem();
        end

    end

    methods (Access = private)

        function init(obj)
            close all
        end

        function solveProblem(obj)
            s.mesh = obj.mesh;
            s.initialPhaseField = obj.initialPhaseField;
            s.materialPhaseField = obj.materialPhaseField;
            s.dissipationPhaseField = obj.dissipationPhaseField;
            s.l0 = obj.l0;
            s.Constant = obj.Constant;
            s.bcVal = obj.bcVal;
            PhaseFieldComputer(s);
        end

        function createMesh(obj)
            %obj.createOneElementMesh();
            %obj.createTwoElementMesh();
            obj.createArbitraryElementMesh(20);
            %obj.createFiberMatrixMesh();
            %obj.createSingleEdgeNotchedMesh();
            %obj.createLshapeMesh();
        end

        function createInitialPhaseField(obj)
            sAF.fHandle = @(x) 0*x;
            sAF.ndimf   = 1;
            sAF.mesh    = obj.mesh;
            xFun = AnalyticalFunction(sAF);

            phi = xFun.project('P1');
            obj.initialPhaseField = phi;
        end

        function createMaterialPhaseField(obj)
            s.mesh = obj.mesh;
            s.materialInterpolation = obj.createMaterialInterpolation();
            s.E = obj.E;
            s.nu = obj.nu;
            s.Gc = obj.Gc;
            %s.fc = obj.fc;
            obj.materialPhaseField = MaterialPhaseField(s);
        end

        function createDissipationInterpolation(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation = 'PhaseFieldD';
            s.pExp = obj.pExp;
            obj.dissipationPhaseField = MaterialInterpolation.create(s);

            if s.pExp == 1
                obj.Constant = obj.Gc/(4*(2/3));
            elseif s.pExp == 2
                obj.Constant = obj.Gc/(4*(1/2));
            end
        end

    end

    methods (Access = private)

        function createOneElementMesh(obj)
            sM.coord = [0,0;
                1,0;
                1,1;
                0,1];
            sM.connec = [1 2 3 4];

            m = Mesh(sM);
            m.plot();
            obj.mesh = m;
        end

        function createTwoElementMesh(obj)
            sM.coord = [0,0;
                        1,0;
                        2,0;
                        0,1;
                        1,1;
                        2,1];
            sM.connec = [1 2 5 4
                2 3 6 5];

            m = Mesh(sM);
            m.plot();
            obj.mesh = m;
        end

        function createArbitraryElementMesh(obj,n)
            m = UnitQuadMesh(n,n/4);
            obj.mesh = m;
        end

        function createFiberMatrixMesh(obj)
            % Generate coordinates
            x1 = linspace(0,1,20);
            x2 = linspace(1,2,20);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            sBg.coord  = V(:,1:2);
            sBg.connec = F;
            bgMesh = Mesh(sBg);
            bdMesh  = bgMesh.createBoundaryMesh();
            
            % Level set creation
            sLS.type       = 'circleInclusion';
            sLS.mesh       = bgMesh;
            sLS.ndim       = 2;
            sLS.fracRadius = 0.4;
            sLS.coord      = bgMesh.coord;
            ls = LevelSetCreator.create(sLS);
            levelSet = ls.getValue();

            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bdMesh;
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
            
            obj.mesh = uMesh.createInnerMeshGoodConditioning();
            obj.mesh.plot;
        end

        function createSingleEdgeNotchedMesh(obj)
            file = 'PF_SENmesh';
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function createLshapeMesh(obj)
            file = 'PF_Lmesh';
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function matInt = createMaterialInterpolation(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation = 'PhaseFieldI';
            s.nElem = obj.mesh.nelem;
            s.dim = '2D';
            s.constitutiveProperties.E = obj.E;
            s.constitutiveProperties.nu = obj.nu;

            matInt = MaterialInterpolation.create(s);
        end
    end

end