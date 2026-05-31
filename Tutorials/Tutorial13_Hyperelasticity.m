classdef Tutorial13_Hyperelasticity < handle
    
    properties (Access = public)
        output
    end

    properties (Access = private)
        mesh
        boundaryConditions
        material, matProp
        functional

        eifemData

        tolSameNode
        nSubdomains

        referenceMesh
        bS
        iC
        lG
        iCR
        discMesh
    end

    methods (Access = public)

        function obj = Tutorial13_Hyperelasticity()
            obj.init()
            obj.createMesh();
            obj.createBoundaryConditions()
            obj.createMaterial()
            obj.createFunctional()

            obj.createEIFEMdata()

            obj.solveHyperelasticityProblem()
        end

        function solveHyperelasticityProblem(obj)

            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.functional         = obj.functional;
            s.material.tensor    = obj.material;
            s.material.prop      = obj.matProp;

            s.monitoring.set       = false;
            s.monitoring.printInfo = true;
            s.monitoring.printFile = false;
            s.monitoring.fileNameOut = 'NeoElastic';
            s.tolerance = 1e-12;
            s.maxIter   = 100;

            % % EIFEM
            % s.eifemData = obj.eifemData;
            % s.activePreconditioner = 'PCG_EIFEM';
            % s.compareEIFEM = true;
            % %

            % EIFEM / ILU constants
            s.eifemData = obj.eifemData;
            
            % Primer pas: Precondicionadors constants
            s.activePreconditioner = 'PCG_ILU_EIFEM_ILU_CONSTANT';
            
           
            
            s.compareEIFEM = true;
            s.useConstantPreconditioners = true;
            %

            hyperComp = HyperelasticityComputer(s);
            obj.output = hyperComp.compute();

        end


    end

    methods (Access = private)

        function init(obj)
            close all

            %
            obj.tolSameNode = 1e-10;
            %obj.nSubdomains = [15 5];
            obj.nSubdomains = [2 1];
        end

        function createMesh(obj)
            % meshType = 'Hole';
            % switch meshType
            %     case {'Hole', 'HoleDirich'}
            %         IM = Mesh.createFromGiD('holeMeshQuad.m');
            %         obj.mesh = IM;
            %     case {'Bending', 'Traction'}
            %         obj.mesh = UnitQuadMesh(20,20);
            %     case {'Metamaterial'}
            %         load('NegativePoissonMesh.mat','NegPoissMesh');
            %         s.coord  = NegPoissMesh.coord;
            %         s.connec = NegPoissMesh.connec;
            %         obj.mesh = Mesh.create(s);
            %     otherwise
            %         obj.mesh = QuadMesh(10,1,100,100);
            %         %obj.mesh = HexaMesh(2,1,1,20,5,5);
            % end

            % Reference mesh d'EIFEM
            load('DEF_Q4porL_1.mat');
            sR.coord     = EIFEoper.MESH.COOR;
            sR.connec    = EIFEoper.MESH.CN;
            sR.interType = 'QUADRATIC';

            % CANVI 1: Correcció nodes cantonada (copiat de TutorialEIFEM.createReferenceMesh)
            tol = 1e-8;
            xmax = max(sR.coord(:,1)); xmin = min(sR.coord(:,1));
            ymax = max(sR.coord(:,2)); ymin = min(sR.coord(:,2));
            
            mask = abs(sR.coord(:,1) - xmax) < tol & abs(sR.coord(:,2) - ymax) < tol;
            sR.coord(mask, :) = sR.coord(mask, :) - [1e-9, 0];
            
            mask = abs(sR.coord(:,1) - xmax) < tol & abs(sR.coord(:,2) - ymin) < tol;
            sR.coord(mask, :) = sR.coord(mask, :) - [1e-9, 0];
            
            mask = abs(sR.coord(:,1) - xmin) < tol & abs(sR.coord(:,2) - ymax) < tol;
            sR.coord(mask, :) = sR.coord(mask, :) + [1e-9, 0];
            
            mask = abs(sR.coord(:,1) - xmin) < tol & abs(sR.coord(:,2) - ymin) < tol;
            sR.coord(mask, :) = sR.coord(mask, :) + [1e-9, 0];
            % FI CANVI 1

            obj.referenceMesh = Mesh.create(sR);
        
            obj.bS = obj.referenceMesh.createBoundaryMesh();
        
            s.nsubdomains   = obj.nSubdomains;
            s.meshReference = obj.referenceMesh;
            s.tolSameNode   = obj.tolSameNode;
        
            m = MeshCreatorFromRVE.create(s);
            [obj.mesh,~,obj.iC,~,obj.lG,obj.iCR,obj.discMesh] = m.create();

        end

        function createBoundaryConditions(obj)
            s.type = 'DisplacementTractionX';
            %s.type = 'ForceTractionXClamped';
            s.values = linspace(0,1,101);
            obj.boundaryConditions = BoundaryConditionsCreator(obj.mesh,s);
        end

        function createMaterial(obj)
            mu     = ConstantFunction.create(1,obj.mesh);
            lambda = ConstantFunction.create(1,obj.mesh);
            ndim   = obj.mesh.ndim;
            kappa = IsotropicElasticMaterial.computeKappaFromShearAndLambda(mu,lambda,ndim);

            s.type  = 'ISOTROPIC';
            s.ndim  = ndim;
            s.bulk  = kappa;
            s.shear = mu;
            obj.material = Material.create(s);
            obj.matProp.mu = mu;
            obj.matProp.lambda = lambda;
        end

        function createFunctional(obj)
            obj.matProp.type = 'Neohookean';
            s.matTensor    = obj.material;
            s.quadOrder    = 3;
            s.matProp      = obj.matProp;
            s.mesh         = obj.mesh;
            s.testSpace.u  = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.functional = ElasticityFunctional(s);
        end

        function createEIFEMdata(obj)

            minx = min(obj.mesh.coord(:,1));
            maxx = max(obj.mesh.coord(:,1));
        
            isLeft  = @(coor) (abs(coor(:,1) - minx) < obj.tolSameNode);
            isRight = @(coor) (abs(coor(:,1) - maxx) < obj.tolSameNode);
        
            Dir{1}.domain    = @(coor) isLeft(coor);
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;
        
            Dir{2}.domain    = @(coor) isRight(coor);
            Dir{2}.direction = [1,2];
            Dir{2}.value     = 0;
        
            dirichletFun = [];

            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.mesh, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end
        
            sBC.pointloadFun = [];
            sBC.dirichletFun = dirichletFun;
            sBC.periodicFun  = [];
            sBC.mesh         = obj.mesh;
            bc               = BoundaryConditions(sBC);
        
            sAp.mesh               = obj.mesh;
            sAp.boundaryConditions = bc;
            bcApplier              = BCApplier(sAp);
            

            obj.eifemData.mesh          = obj.mesh;
            obj.eifemData.fileNameEIFEM = 'DEF_Q4porL_1.mat';
            obj.eifemData.referenceMesh = obj.referenceMesh;
            obj.eifemData.dir           = Dir;
            obj.eifemData.iC            = obj.iC;
            obj.eifemData.lG            = obj.lG;
            obj.eifemData.bS            = obj.bS;
            obj.eifemData.iCR           = obj.iCR;
            obj.eifemData.discMesh      = obj.discMesh;
            obj.eifemData.nSubdomains   = obj.nSubdomains;
            obj.eifemData.tolSameNode   = obj.tolSameNode;
            obj.eifemData.bcApplier     = bcApplier;
            obj.eifemData.nNodes        = obj.mesh.nnodes;
            obj.eifemData.nDimf         = obj.mesh.ndim;
        end

    end

end