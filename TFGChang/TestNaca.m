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
        % AOAd
        % xCentral
        % yCentral
    end
    
    properties (Access = private)
        refMesh
        levelSet
        rawMesh
        mesh
        uMesh
        material
        velocityFun
        pressureFun
        forcesFormula
        dirConditions
        dirDofs
        nodesConditions
    end
    
    methods (Access = public)
        
        function obj = TestNaca(cParams)
            close all;
            obj.init(cParams);
            obj.createReferenceMesh();
            obj.verifyReferenceMesh();
            obj.levelSet = obj.createLevelSet(obj.refMesh, obj.p, obj.M, obj.t);
            obj.createFluidMesh();
            obj.createFluidMeshGoodConditioning();
            obj.verifyFinalMesh();
            obj.createMaterial();
            obj.checkMaterial();
            obj.createTrialFunction();
            obj.defineBoundaryConditions();
            obj.checkBoundaryConditions();
        end

        function compute(obj)
            obj.solveStokesProblem();
            obj.verifyStokeProblemSolver();
            obj.plotResults();
            obj.CalculateAeroForces();
            obj.verifyAeroForces();
        end

        function validate(obj)
            obj.verifyReferenceMesh();
            obj.verifyFinalMesh();
            obj.checkBoundaryConditions();
            obj.verifyStokeProblemSolver();
            obj.verifyAeroForces();
        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.p        = cParams.p;
            obj.M        = cParams.M;
            obj.t        = cParams.t;
            % obj.AOAd     = cParams.AOAd;
            % obj.xCentral = cParams.xCentral;
            % obj.yCentral = cParams.yCentral;
        end
                
        function createReferenceMesh(obj)      
            obj.refMesh = TriangleMesh(2,1,150,75);
            %QuadMesh(10,4,100,40);
        end
        
        function createFluidMesh(obj)
            s.backgroundMesh = obj.refMesh;
            s.boundaryMesh   = obj.refMesh.createBoundaryMesh();
            obj.uMesh            = UnfittedMesh(s);
            obj.uMesh.compute(obj.levelSet.fValues);       
            obj.rawMesh = obj.uMesh.createInnerMesh();
        end

        function createFluidMeshGoodConditioning(obj)
            m        = obj.createMeshAlphaTriangulation();
            s.connec = obj.computeConnectivitiesGoodCond(m);
            s.coord  = obj.rawMesh.coord;
            m2       = Mesh.create(s);
            obj.mesh = m2.computeCanonicalMesh();
            figure;
            obj.mesh.plot();
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
            xV(:,1,:) = mAlpha.computeBaricenter();
            lsFun  = obj.createLevelSet(mAlpha, obj.p, obj.M, obj.t);
            lsElem = squeeze(lsFun.evaluate(xV));
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
            s.mesh         = obj.mesh;
            s.uMesh        = obj.uMesh;
            s.velocityFun  = obj.velocityFun;
            s.pressureFun  = obj.pressureFun;
            BCClass              = BoundaryCondition(s); 
            BCClass.compute();
            obj.forcesFormula    = BCClass.forcesFormula;
            obj.dirConditions    = BCClass.dirConditions;
            obj.dirDofs          = BCClass.dirDofs;
            obj.nodesConditions  = BCClass.nodesConditions;
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
            obj.pressureFun.plot();
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

        function verifyReferenceMesh(obj)
            s.refMesh         = obj.refMesh;
            referenceMeshTest = ReferenceMeshTestComputer(s);
            referenceMeshTest.compute();
        end

        function verifyFinalMesh(obj)
            s.mesh         = obj.mesh;
            finalMeshTest = FinalMeshTestComputer(s);
            finalMeshTest.compute();
        end

        function checkMaterial(obj)
            s.material   = obj.material;
            MaterialTest = MaterialTestComputer(s);
            MaterialTest.compute();
        end

        function checkBoundaryConditions(obj)
            s.forcesFormula      = obj.forcesFormula;
            s.dirConditions      = obj.dirConditions;
            s.dirDofs            = obj.dirDofs;
            s.nodesConditions    = obj.nodesConditions;
            BoundaryConditionsTest = BoundaryConditionsTestComputer(s);
            BoundaryConditionsTest.compute();
        end

        function verifyStokeProblemSolver(obj)
            s.velocityFun       = obj.velocityFun;
            s.pressureFun       = obj.pressureFun;
            StokesProblemSolverTest = StokesProblemSolverTestComputer(s);
            StokesProblemSolverTest.compute();
        end

        function verifyAeroForces(obj)
            s.L    = obj.L;
            s.D    = obj.D;
            AeroForcesCalculationTest = AeroForcesCalculationTestComputer(s);
            AeroForcesCalculationTest.compute();
        end


    end

    methods (Static, Access = private)

        function ls = createLevelSet(m, p, M, t)
            s.type = 'Naca';
            s.xLE  = 0.5;
            s.yLE  = 0.5;

            s.chord = 1;
            s.p     = p;
            s.m     = M;
            s.t     = t;

            g  = GeometricalFunction(s);
            ls = g.computeLevelSetFunction(m);
        end

    end

end