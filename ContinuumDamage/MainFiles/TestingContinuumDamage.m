classdef TestingContinuumDamage < handle

    properties (Access = public)
        data
    end

    properties (Access = private)
        benchmark
        tolerance
    end

    properties (Access = private)
        mesh
        boundaryConditions
        internalDamageVariable
        damageFunctional
    end

    methods (Access = public)

        function obj = TestingContinuumDamage()
            obj.init();
            obj.mesh                   = obj.createMesh();
            obj.boundaryConditions     = obj.createBoundaryConditions();
            obj.damageFunctional       = obj.createContinuumDamageFunctional();
            obj.internalDamageVariable = obj.createInternalDamageVariable();
        end

        function outputData = compute(obj)
            s.mesh                   = obj.mesh;
            s.boundaryConditions     = obj.boundaryConditions;
            s.damageFunctional       = obj.damageFunctional;
            s.internalDamageVariable = obj.internalDamageVariable;
            s.tol                    = obj.tolerance;
            CDComp = ContinuumDamageComputer(s);

            outputData = CDComp.compute();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.benchmark = cParams.benchmark;
        end

        function mesh = createMesh(obj)
            if obj.benchmark.mesh.type == "Rectangle"
                l = obj.benchmark.mesh.length;
                w = obj.benchmark.mesh.width;
                N = obj.benchmark.mesh.lN;
                M = obj.benchmark.mesh.wN;
                mesh = QuadMesh(l,w,N,M);
            else
                file = obj.benchmark.mesh.type;
                a.fileName = file;
                s = FemDataContainer(a);
                mesh = s.mesh;
            end
        end

        function bc = createBoundaryConditions(obj)
            s.mesh = obj.mesh;
            s.bcType = obj.benchmark.bc.type;
            s.bcValueSet = obj.benchmark.bc.bcValues;            
            bc = BcContinuumDamage(s); %% Merge w/ PFBoundaryCreator -> BenchmarkBoundaryCreator
        end

        function r = createInternalDamageVariable(obj)
            s.r0 = obj.createR0();
            s.mesh = obj.mesh;
            r = InternalDamageVariable(s);
        end

        function sF = createContinuumDamageFunctional(obj)
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.material           = obj.createDamagedMaterial();
            s.quadOrder          = 2;
            s.test               = LagrangianFunction.create(obj.mesh,2,'P1');
            sF = ShFunc_ContinuumDamage(s);
        end        
    
        function dM = createDamagedMaterial(obj)
            s.baseMaterial = obj.createBaseMaterial();
            s.damage       = obj.createDamagedLaw();
            dM = DamagedMaterial(s);
        end

        function mat = createBaseMaterial(obj)
            E = 210;
            nu = 0.3;
            s.type    = 'ISOTROPIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = ConstantFunction.create(E,obj.mesh);
            s.poisson = ConstantFunction.create(nu,obj.mesh);
            mat = Material.create(s);
        end        

        function d = createDamagedLaw(obj)
            s.hardeningLaw = obj.createHardeningLaw();
            d = DamageLaw(s);
        end

        function hL = createHardeningLaw(obj)
            r1     = 20;
            s.r1   = ConstantFunction.create(r1,obj.mesh);
            s.type = 'Linear';            
            s.H    = -0.5;
            s.A    = 0.1;
            s.r0   = obj.createR0();
            s.w1   = 8;
            s.qInf = 2;
            hL = HardeningLaw.create(s);
        end          

        function r0 = createR0(obj)
            r0 = 10;
            r0 = ConstantFunction.create(r0,obj.mesh);            
        end

 
    end
end