classdef TestingHyperelasticity < handle

    properties (Access = private)
        meshType
        bcProp
        matProp
        monitoring
        tolerance
        maxIter
    end

    properties (Access = private)
        mesh
        boundaryConditions
        material
        functional
    end

    methods (Access = public)

        function obj = TestingHyperelasticity(cParams)
            obj.init(cParams)
            obj.createMesh();
            obj.createBoundaryConditions()
            obj.createMaterial()
            obj.createFunctional()
        end

        function outputData = compute(obj)
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.functional         = obj.functional;
            s.material.tensor    = obj.material;
            s.material.prop      = obj.matProp;
            s.monitoring         = obj.monitoring;
            s.tolerance          = obj.tolerance;
            s.maxIter            = obj.maxIter;
            hyperComp = HyperelasticityComputer(s);
            
            outputData = hyperComp.compute();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.meshType   = cParams.meshType;
            obj.bcProp     = cParams.bcProp;
            obj.matProp    = cParams.matProp;
            obj.monitoring = cParams.monitoring;
            obj.tolerance  = cParams.tolerance;
            obj.maxIter    = cParams.maxIter;
        end

        function createMesh(obj)
            switch obj.meshType
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
            nSteps = obj.bcProp.nSteps;
            maxVal = obj.bcProp.maxVal;
            type   = obj.bcProp.type;
            cParams.values = linspace(0,maxVal,nSteps);
            cParams.type = type;
            bc = BoundaryCreator(obj.mesh,cParams);
            obj.boundaryConditions = bc;
        end

        function createMaterial(obj)
            mu     = ConstantFunction.create(obj.matProp.mu,obj.mesh);
            lambda = ConstantFunction.create(obj.matProp.lambda,obj.mesh);
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
            s.matTensor    = obj.material;
            s.quadOrder    = 3;
            s.matProp      = obj.matProp;
            s.mesh         = obj.mesh;
            s.testSpace.u  = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.functional = ElasticityFunctional(s);
        end

    end

end