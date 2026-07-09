classdef Tutorial13_Hyperelasticity < handle
    
    properties (Access = public)
        output
    end

    properties (Access = private)
        mesh
        boundaryConditions
        material, matProp
        functional
    end

    methods (Access = public)

        function obj = Tutorial13_Hyperelasticity()
            obj.init()
            obj.createMesh();
            obj.createBoundaryConditions()
            obj.createMaterial()
            obj.createFunctional()
            obj.solveHyperelasticityProblem()
        end

        function solveHyperelasticityProblem(obj)
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.functional         = obj.functional;
            s.material.tensor    = obj.material;
            s.material.prop      = obj.matProp;

            s.monitoring.set       = true;
            s.monitoring.printInfo = true;
            s.monitoring.printFile = false;
            s.monitoring.fileNameOut = 'NeoElastic';
            s.tolerance = 1e-12;
            s.maxIter   = 100;

            hyperComp = HyperelasticityComputer(s);
            obj.output = hyperComp.compute();
        end
    end

    methods (Access = private)

        function init(obj)
            close all
        end

        function createMesh(obj)
            meshType = 'Hole';
            switch meshType
                case {'Hole', 'HoleDirich'}
                    IM = Mesh.createFromGiD('holeMeshQuad.m');
                    obj.mesh = IM;
                case {'Bending', 'Traction'}
                    obj.mesh = UnitQuadMesh(20,20);
                case {'Metamaterial'}
                    load('NegativePoissonMesh.mat','NegPoissMesh');
                    s.coord  = NegPoissMesh.coord;
                    s.connec = NegPoissMesh.connec;
                    obj.mesh = Mesh.create(s);
                otherwise
                    obj.mesh = QuadMesh(10,1,100,100);
                    %obj.mesh = HexaMesh(2,1,1,20,5,5);
            end
        end

        function createBoundaryConditions(obj)
            s.type = 'DisplacementTractionX';
            s.values = linspace(0,1,21);
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

    end

end