classdef TestNaca < handle
    
    properties (Access = public)
        L
        D
        E
        velocityFun
        pressureFun
    end
    
    properties (Access = private)
        M
        p
        t
        AoA
        chord
        length
        height
        nx
        ny
        flowType
        uRef
    end
    
    properties (Access = private)
        refMesh
        levelSet
        rawMesh
        mesh
        backupmesh
        uMesh
        material
        filter
        forcesFormula
        dirConditions
        dirDofs
        nodesConditions
    end
    
    methods (Access = public)
        
        function obj = TestNaca(cParams)
            %close all;
            obj.init(cParams);    
            obj.createReferenceMesh();
            [AirfoilParams, BGParams] = obj.setParams();
            obj.createLevelSet(AirfoilParams, BGParams);
            obj.createFluidMesh();
            %obj.plotMesh();
            %obj.createFluidMeshGoodConditioning(AirfoilParams, BGParams));
            obj.createMaterial();
            obj.createTrialFunction();
            obj.defineBoundaryConditions();
            obj.defineAppliedForces();
        end

        function compute(obj)
            obj.createPressureFilter();
            obj.solveProblem();
            obj.plotResults();
            obj.CalculateAeroForces();
        end

        function validate(obj)
            obj.verifyAeroForces();
        end

        function print(obj)
            obj.printResults();
        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.p        = cParams.p;
            obj.M        = cParams.M;
            obj.t        = cParams.t;
            obj.AoA      = cParams.AoA;
            obj.chord    = cParams.chord;
            obj.length   = cParams.length;
            obj.height   = cParams.height;
            obj.nx       = cParams.nx;
            obj.flowType = cParams.flowType;
            obj.uRef     = cParams.uRef;
        end

        function [AirfoilParams, BGParams] = setParams(obj)
            AirfoilParams.p      = obj.p;
            AirfoilParams.M      = obj.M;
            AirfoilParams.t      = obj.t;
            AirfoilParams.AoA    = obj.AoA;
            AirfoilParams.chord  = obj.chord;
            BGParams.length      = obj.length;
            BGParams.height      = obj.height;       
            BGParams.nx          = obj.nx;
        end
                
        function createReferenceMesh(obj)
            % obj.length  = 8;
            % obj.height  = 4;
            % nx          = 300;
            % ny          = 150;
            % obj.refMesh = QuadMesh(obj.length,obj.height,nx,ny); 
            obj.ny      = round(obj.nx / obj.length * obj.height / 0.8);
            obj.refMesh = TriangleMesh(obj.length,obj.height,obj.nx,obj.ny);
        end

        function createLevelSet(obj,AirfoilParams, BGParams)
            g = obj.createNacaFunction(AirfoilParams, BGParams);
            obj.levelSet = g.computeLevelSetFunction(obj.refMesh);
            %obj.levelSet.plot();
        end
   
        function createFluidMesh(obj)
            s.backgroundMesh = obj.refMesh;
            s.boundaryMesh   = obj.refMesh.createBoundaryMesh();
            obj.uMesh        = UnfittedMesh(s);
            obj.uMesh.compute(obj.levelSet.fValues);       
            obj.rawMesh = obj.uMesh.createInnerMesh();
            obj.mesh    = obj.rawMesh;
        end

        function plotMesh(obj)
            figure;
            obj.mesh.plot();
            title("Mesh with the airfoil inclusion.");
            xlabel("x");
            ylabel("y");
        end

        function createFluidMeshGoodConditioning(obj,AirfoilParams, BGParams)
            m              = obj.createMeshAlphaTriangulation();
            s.connec       = obj.computeConnectivitiesGoodCond(m,AirfoilParams, BGParams);
            s.coord        = obj.rawMesh.coord;
            m2             = Mesh.create(s);
            obj.backupmesh = m2.computeCanonicalMesh();
            %figure;
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

        function connec = computeConnectivitiesGoodCond(obj,mAlpha,AirfoilParams, BGParams)
            q        = Quadrature.create(mAlpha, 0);
            xV       = q.posgp;
            g  = obj.createNacaFunction(AirfoilParams, BGParams);
            lsElem = squeeze(g.evaluate(xV,mAlpha));
            connec = mAlpha.connec(lsElem<=0,:);
        end

        function createMaterial(obj)
            e.type       = obj.flowType;  
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
            s.uRef         = obj.uRef;
            BCClass              = StokesProblemBoundaryCondition(s); 
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

        function createPressureFilter(obj)
            s.filterType = 'PDE';
            s.mesh       = obj.mesh;
            s.trial      = obj.pressureFun;
            obj.filter   = Filter.create(s);
        end

        function s = initProblemSolver(obj)
            s.mesh             = obj.mesh;
            s.forcesFormula    = obj.forcesFormula;
            s.dirConditions    = obj.dirConditions;
            s.dirDofs          = obj.dirDofs;
            s.material         = obj.material;
            s.velocityFun      = obj.velocityFun;
            s.pressureFun      = obj.pressureFun;
        end

        function solveProblem(obj)
            s = obj.initProblemSolver();
        
            switch obj.flowType
                case 'Stokes'
                    obj.solveStokesProblem(s);
                case 'NavierStokes'
                    s.velocityField = LagrangianFunction.create(obj.mesh, 2, 'P2');
                    obj.solveNavierStokesProblem(s);
                otherwise
                    error('Unknown flow type: %s', obj.flowType);
            end
           
            obj.pressureFun = obj.filter.compute(obj.pressureFun, 3);

        end

        function solveStokesProblem(obj,s)
            SolverResults      = StokesProblemSolver(s); 
            SolverResults.compute();
            obj.velocityFun    = SolverResults.velocityFun;
            obj.pressureFun    = SolverResults.pressureFun;
        end

        function solveNavierStokesProblem(obj,s)
            s.velocityField    = LagrangianFunction.create(obj.mesh, 2, 'P2');
            SolverResults      = NavierStokesProblemSolver(s); 
            SolverResults.compute();
            obj.velocityFun    = SolverResults.velocityFun;
            obj.pressureFun    = SolverResults.pressureFun;
        end

        function plotResults(obj)     
            TestNaca.plotVelocity(obj.velocityFun);
            TestNaca.plotPressure(obj.pressureFun);
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

        function printResults(obj)
            fileID = fopen('results.txt', 'a');
            fprintf(fileID, '%.2f %.2f %.2f %.2f %.4f %.4f %.4f\n', obj.M, obj.p, obj.t, obj.AoA, obj.L, obj.D, obj.E);
            fclose(fileID);
        end


    end

    methods (Static, Access = public)

        function plotVelocity(velocityFun)
            v = LagrangianFunction.create(velocityFun.mesh, 1, 'P2');
            v.setFValues(sqrt(velocityFun.fValues(:,1).^2 + velocityFun.fValues(:,2).^2));
            v.plot()
            title("Velocity distribution.");
            xlabel("x");
            ylabel("y");
            colorbar;
            axis image tight;
            %set(gca, 'Position', [0.1 0.1 0.75 1]);
            set(gcf, 'WindowState', 'maximized');
        end

        function plotVelocityX(velocityFun)
            vx = LagrangianFunction.create(velocityFun.mesh, 1, 'P2');
            vx.setFValues(velocityFun.fValues(:,1));
            vx.plot();
            title("Velocity distribution in the x direction.");
            xlabel("x");
            ylabel("y");
            colorbar;
            axis equal tight;  
            %set(gca, 'Position', [0.1 0 0.75 1]);
            set(gcf, 'WindowState', 'maximized');
        end

        function plotVelocityY(velocityFun)
            vy = LagrangianFunction.create(velocityFun.mesh, 1, 'P2');
            vy.setFValues(velocityFun.fValues(:,2));
            vy.plot();
            title("Velocity distribution in the y direction."); 
            xlabel("x");
            ylabel("y"); 
            colorbar;
            axis equal tight; 
            %set(gca, 'Position', [0.1 0 0.75 1]);
            set(gcf, 'WindowState', 'maximized');
        end

        function plotPressure(PressureFun)
            PressureFun.plot();
            xlabel("x");
            ylabel("y");
            title("Pressure distribution"); 
            axis equal tight;
            colorbar;
            caxis([-20, 20]);
            %set(gca, 'Position', [0.1 0 0.75 1]);
            set(gcf, 'WindowState', 'maximized');
        end

    end

    methods (Static, Access = private)

        function g = createNacaFunction(AirfoilParams, BGParams)
            s.type = 'Naca'; % 'Naca
            s.xLE  = (BGParams.length - AirfoilParams.chord) / 2;
            s.yLE  = BGParams.height/2;

            s.chord = AirfoilParams.chord;
            s.p     = AirfoilParams.p;
            s.m     = AirfoilParams.M;
            s.t     = AirfoilParams.t;
            s.AoA   = AirfoilParams.AoA;

            g  = GeometricalFunction(s);
        end

    end

    

end