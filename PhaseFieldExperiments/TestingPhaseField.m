classdef TestingPhaseField < handle

    properties (Access = public)
        E  = 210;
        nu = 0.3;
        Gc = 5e-3;
        l0 = 0.1;
        pExp = 2;
        % bcVal = [linspace(0,8.5e-3,100), ...
        %         linspace(8.5e-3,0,100), ...
        %         linspace(0,-1e-2,100), ...
        %         ];
        bcVal = linspace(1e-4,1e-1,100);
        % bcVal = 1;
        % bcVal = [0.001];

        % 1) CHANGE MESH
        % 2) CHANGE BC (TYPE AND DIRECTION)
        % 3) CHANGE REACTIONS (TYPE AND DIRECTION)
    end

    properties (Access = private)
    end

    properties (Access = private)
        mesh
        initialPhaseField
        materialPhaseField
        dissipationPhaseField
        constant
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
            s.constant = obj.constant;
            s.bcVal = obj.bcVal;
            PhaseFieldComputer(s);
        end

        function createMesh(obj)
            %obj.createOneElementMesh();
            %obj.createTwoElementMesh();
            obj.createArbitraryElementMesh(1,1,10,10);
            %obj.createFiberMatrixMesh();
            %obj.createSingleEdgeNotchedMesh();
            %obj.createLshapeMesh();
        end

        function createInitialPhaseField(obj)
            phi = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.initialPhaseField = phi;
        end

        function createMaterialPhaseField(obj)
             sParam.mesh = obj.mesh;
             sParam.order = 'P1';
             sParam.fValues = ones(obj.mesh.nnodes,1)*obj.E;
             s.young = LagrangianFunction(sParam);
             sParam.fValues = ones(obj.mesh.nnodes,1)*obj.nu;
             s.poisson = LagrangianFunction(sParam);
             
             sIso.ndim = obj.mesh.ndim;
             sIso.young = s.young ;
             sIso.poisson = s.poisson;

             
             s.mesh         = obj.mesh;
             s.baseMaterial = Isotropic2dElasticMaterial(sIso);  
             s.degradation  = obj.createMaterialInterpolation();
             s.Gc = obj.Gc;
             
             obj.materialPhaseField = MaterialPhaseField(s);

            % sH.fileName    = 'IsoMicroDamage';
            % obj.materialPhaseField = HomogenizedPhaseField(sH);
        end

        function createDissipationInterpolation(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation = 'PhaseFieldD';
            s.pExp = obj.pExp;
            obj.dissipationPhaseField = MaterialInterpolator.create(s);

            if s.pExp == 1
                obj.constant = obj.Gc/(4*(2/3));
            else%if s.pExp == 2
                obj.constant = obj.Gc/(4*(1/2));
            end
        end

        function matInt = createMaterialInterpolation(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation = 'PhaseFieldI';
            matInt = MaterialInterpolator.create(s);
        end

    end

    methods (Access = private)

        function createOneElementMesh(obj)
            sM.coord = [0,0;
                1,0;
                1,1;
                0,1];
            sM.connec = [1 2 3 4];

            m = Mesh.create(sM);
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

            m = Mesh.create(sM);
            m.plot();
            obj.mesh = m;
        end

        function createArbitraryElementMesh(obj,lx, ly, nx,ny)
            m = QuadMesh(lx, ly, nx, ny);
            m.plot();
            obj.mesh = m;
        end

        function createFiberMatrixMesh(obj) %%%%%% REVIEW LEVEL SETS %%%%%%%%%%%%%%%
            % Generate coordinates
            x1 = linspace(0,1,20);
            x2 = linspace(1,2,20);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            sBg.coord  = V(:,1:2);
            sBg.connec = F;
            bgMesh = Mesh.create(sBg);
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
    end

end