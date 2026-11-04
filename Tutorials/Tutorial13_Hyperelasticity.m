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
            s.tolerance = 1e-8;
            s.maxIter   = 100;
            hyperComp = HyperelasticityComputer(s);
            obj.output = hyperComp.compute();
        end
    end
    methods (Access = private)
        function init(obj)
            close all
        end
        function mesh = createRefMesh(obj)
            meshType = 'Metamaterial';
            switch meshType
                case {'Hole', 'HoleDirich'}
                    IM = Mesh.createFromGiD('holeMeshQuad.m');
                    mesh = IM;
                case {'Bending', 'Traction'}
                    mesh = UnitQuadMesh(20,20);
                case {'Metamaterial'}
                    load('NegativePoissonMesh.mat','NegPoissMesh');
                    s.coord  = NegPoissMesh.coord;
                    s.connec = NegPoissMesh.connec;
                    mesh = Mesh.create(s);
                otherwise
                    mesh = QuadMesh(10,1,100,100);
                    %obj.mesh = HexaMesh(2,1,1,20,5,5);
            end
        end

        function createMesh(obj)
            refMesh         = obj.createRefMesh();
            s.nsubdomains   = [4,4];
            s.meshReference = refMesh;
            s.tolSameNode   = 1e-8;
            m = MeshCreatorFromRVE.create(s);
            [mD,~,~,~,~,~,~] = m.create();
            obj.mesh = mD;
        end
        function createBoundaryConditions(obj)
            s.type = 'ForceTractionXClamped';
            s.values = linspace(2,4,300);
            obj.boundaryConditions = BoundaryConditionsCreator(obj.mesh,s);
        end
        function createMaterial(obj)
            mu     = ConstantFunction.create(1,obj.mesh);
            lambda = ConstantFunction.create(1,obj.mesh);
            ndim   = obj.mesh.ndim;

            muE     = Expand(mu,4);
            lambdaE = Expand(lambda,4);
            I       = ConstantFunction.create(eye4D(ndim),obj.mesh);
            IxI     = ConstantFunction.create(kronEye(ndim),obj.mesh);
            obj.material = 2*muE.*I + lambdaE.*IxI;

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