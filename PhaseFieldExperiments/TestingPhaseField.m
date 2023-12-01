classdef TestingPhaseField < handle

    properties (Access = public)
        E = 210;
        nu = 0.3;
        Gc = 5e-3;
        fc = 1;
    end

    properties (Access = private)

    end

    properties (Access = private)
        mesh
        initialPhaseField
        materialPhaseField
        dissipationPhaseField
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
            PhaseFieldComputer(s);
        end

        function createMesh(obj)
            obj.createOneElementMesh();
            %obj.createTwoElementMesh();
            %obj.createArbitraryElementMesh(2);
            %obj.createOpenHoleMesh();
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
            s.fc = obj.fc;
            obj.materialPhaseField = MaterialPhaseField(s);
        end

        function createDissipationInterpolation(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation = 'PhaseFieldD';
            obj.dissipationPhaseField = MaterialInterpolation.create(s);
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
            sM.coord = [-1,-1;
                1,-1;
                1,1;
                -1,1;
                -1 3;
                1 3];
            sM.connec = [1 2 3 4
                4 3 6 5];

            m = Mesh(sM);
            m.plot();
            obj.mesh = m;
        end

        function createArbitraryElementMesh(obj,n)
            m = UnitQuadMesh(n,n);
            obj.mesh = m;
        end

        function createOpenHoleMesh(obj)
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