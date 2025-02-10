classdef TestNaca < handle
    
    properties (Access = public)
        L
        D
        E
    end
    
    properties (Access = private)
        M
        p
        t
        AoA
        % xCentral
        % yCentral
    end
    
    properties (Access = private)
        refMesh
        levelSet
        rawMesh
        mesh
        backupmesh
        uMesh
        material
        velocityFun
        pressureFun
        forcesFormula
        dirConditions
        dirDofs
        nodesConditions
        length
        height
    end
    
    methods (Access = public)
        
        function obj = TestNaca(cParams)
            close all;
            obj.init(cParams);
            obj.createReferenceMesh();
            obj.createLevelSet();
            obj.createFluidMesh();
            %obj.createFluidMeshGoodConditioning();
            obj.createMaterial();
            obj.createTrialFunction();
            obj.defineBoundaryConditions();
            obj.defineAppliedForces();
        end

        function compute(obj)
            obj.solveStokesProblem();
            obj.plotResults();
            obj.CalculateAeroForces();
        end

        function validate(obj)
            obj.verifyAeroForces();
        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.p        = cParams.p;
            obj.M        = cParams.M;
            obj.t        = cParams.t;
            obj.AoA      = cParams.AoA;
            % obj.xCentral = cParams.xCentral;
            % obj.yCentral = cParams.yCentral;
        end
                
        function createReferenceMesh(obj)
            % obj.length  = 10;
            % obj.height  = 4;
            % nx          = 100;
            % ny          = 40;
            % obj.refMesh = QuadMesh(obj.length,obj.height,nx,ny); 
            obj.length  = 2;
            obj.height  = 1;
            nx          = 150;
            ny          = 75;
            obj.refMesh = TriangleMesh(obj.length,obj.height,nx,ny);
        end

        function createLevelSet(obj)
            g = obj.createNacaFunction(obj.p, obj.M, obj.t, obj.AoA);
            obj.levelSet = g.computeLevelSetFunction(obj.refMesh);
        end
        
        function createFluidMesh(obj)
            s.backgroundMesh = obj.refMesh;
            s.boundaryMesh   = obj.refMesh.createBoundaryMesh();
            obj.uMesh        = UnfittedMesh(s);
            obj.uMesh.compute(obj.levelSet.fValues);       
            obj.rawMesh = obj.uMesh.createInnerMesh();
            obj.mesh    = obj.rawMesh;
            %obj.mesh.plot();
            title("Mesh with the airfoil inclusion.");
            xlabel("x");
            ylabel("y");
            obj.uMesh.plot();
        end

        function createFluidMeshGoodConditioning(obj)
            m              = obj.createMeshAlphaTriangulation();
            s.connec       = obj.computeConnectivitiesGoodCond(m);
            s.coord        = obj.rawMesh.coord;
            m2             = Mesh.create(s);
            obj.backupmesh = m2.computeCanonicalMesh();
            figure;
            %obj.backupmesh.plot();
        end

        function m = createMeshAlphaTriangulation(obj)
            points = obj.rawMesh.coord;
            r      = obj.computeAlphaDistance();
            T      = alphaShape(points,r);
            DT     = alphaTriangulation(T);

            s.connec = DT;
            s.coord  = points;
            m        = Mesh.create(s);
        end

        function r = computeAlphaDistance(obj)
            bCMesh = obj.uMesh.boundaryCutMesh.mesh;
            x      = bCMesh.coord(:,1);
            y      = bCMesh.coord(:,2);
            xMax   = max(x);
            xMin   = min(x);
            yMax   = max(y);
            yMin   = min(y);
            r      = 0.25*min(xMax-xMin,yMax-yMin);
            h      = obj.refMesh.computeMeanCellSize();
            if h>=r
                warning('Alpha radius is too small');
            end
        end

        function connec = computeConnectivitiesGoodCond(obj,mAlpha)
            q  = Quadrature.create(mAlpha, 0);
            xV = q.posgp;
            g  = obj.createNacaFunction(obj.p, obj.M, obj.t, obj.AoA);
            lsElem = squeeze(g.evaluate(xV,mAlpha));
            connec = mAlpha.connec(lsElem<=0,:);
        end

        function createMaterial(obj)
            e.type       = 'STOKES'; 
            e.nelem      = obj.mesh.nelem;
            obj.material = Material.create(e); 
        end

        function createTrialFunction(obj)
            obj.velocityFun = LagrangianFunction.create(obj.mesh, 2, 'P2');
            obj.pressureFun = LagrangianFunction.create(obj.mesh, 1, 'P1');
        end

        function defineBoundaryConditions(obj)
            s.height       = obj.height;
            s.mesh         = obj.mesh;
            s.uMesh        = obj.uMesh;
            s.velocityFun  = obj.velocityFun;
            s.pressureFun  = obj.pressureFun;
            BCClass              = BoundaryCondition(s); 
            BCClass.compute();
            obj.dirConditions    = BCClass.dirConditions;
            obj.dirDofs          = BCClass.dirDofs;
            obj.nodesConditions  = BCClass.nodesConditions;
        end

        function defineAppliedForces(obj)
            sAF.fHandle       = @(coor) [0.*coor,0.*coor];
            sAF.ndimf         = 2;
            sAF.mesh          = obj.mesh;
            obj.forcesFormula = AnalyticalFunction(sAF);
        end

        function solveStokesProblem(obj)
            s.mesh             = obj.mesh;
            s.forcesFormula    = obj.forcesFormula;
            s.dirConditions    = obj.dirConditions;
            s.dirDofs          = obj.dirDofs;
            s.nodesConditions  = obj.nodesConditions;
            s.material         = obj.material;
            s.velocityFun      = obj.velocityFun;
            s.pressureFun      = obj.pressureFun;
            SolverResults      = StokesProblemSolver(s); 
            SolverResults.compute();
            obj.velocityFun    = SolverResults.velocityFun;
            obj.pressureFun    = SolverResults.pressureFun;
        end

        function plotResults(obj)     
            obj.velocityFun.plot(); 
            ax = findall(groot, 'Type', 'axes');
            xlabel(ax(2),"x");
            ylabel(ax(2),"y");
            xlabel(ax(1),"x");
            ylabel(ax(1),"y"); 
            title(ax(2), "Velocity distribution in the x direction.");
            title(ax(1), "Velocity distribution in the y direction.");        
            obj.pressureFun.plot();
            xlabel("x");
            ylabel("y");
            title("Pressure distribution"); 
            caxis([-50 50]);
        end

        function CalculateAeroForces(obj)
            s.mesh             = obj.mesh;
            s.nodesConditions  = obj.nodesConditions;
            s.pressureFun      = obj.pressureFun;
            AeroForcesResults  = AeroForcesCalculation(s); 
            AeroForcesResults.compute();
            obj.L              = AeroForcesResults.L;
            obj.D              = AeroForcesResults.D;
            obj.E              = AeroForcesResults.E;
        end

        function verifyAeroForces(obj)
            s.L    = obj.L;
            s.D    = obj.D;
            AeroForcesCalculationTest = AeroForcesCalculationTestComputer(s);
            AeroForcesCalculationTest.compute();
        end


    end

    methods (Static, Access = private)

        function g = createNacaFunction(p, M, t, AoA)
            s.type = 'Naca';
            s.xLE  = 0.5;
            s.yLE  = 0.5;

            s.chord = 1;
            s.p     = p;
            s.m     = M;
            s.t     = t;
            s.AoA   = AoA;

            g  = GeometricalFunction(s);
        end

    end

end