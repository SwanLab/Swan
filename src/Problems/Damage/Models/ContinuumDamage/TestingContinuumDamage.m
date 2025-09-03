classdef TestingContinuumDamage < handle

    properties (Access = public)
        data
    end

    properties (Access = private)
        benchmark
        matInfo
        damageInfo
        monitoring
        tolerance
        maxIter
    end

    properties (Access = private)
        mesh
        boundaryConditions
        internalDamageVariable
        functional
    end

    methods (Access = public)

        function obj = TestingContinuumDamage(cParams)
            obj.init(cParams);
            obj.mesh                   = obj.createMesh();
            obj.boundaryConditions     = obj.createBoundaryConditions();
            obj.internalDamageVariable = obj.createInternalDamageVariable();
            obj.functional             = obj.createContinuumDamageFunctional();
            
        end

        function outputData = compute(obj)
            s.mesh                   = obj.mesh;
            s.boundaryConditions     = obj.boundaryConditions;
            s.internalDamageVariable = obj.internalDamageVariable;
            s.functional             = obj.functional;
            s.tolerance              = obj.tolerance;
            s.maxIter                = obj.maxIter;
            s.monitoring             = obj.monitoring;
            CDComp = ContinuumDamageComputer(s);

            outputData = CDComp.compute();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.benchmark  = cParams.benchmark;
            obj.matInfo    = cParams.matInfo;
            obj.damageInfo = cParams.damageInfo;
            obj.monitoring = cParams.monitoring;
            obj.tolerance  = cParams.tolerance;
            obj.maxIter    = cParams.maxIter;
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
            s.type     = obj.benchmark.bc.type;
            s.bcValues = obj.benchmark.bc.bcValues;            
            bc = BoundaryCreator(obj.mesh,s);
        end

        function r = createInternalDamageVariable(obj)
            r0     = obj.damageInfo.r0;
            s.r0   = ConstantFunction.create(r0,obj.mesh);
            s.mesh = obj.mesh;
            r = InternalDamageVariable(s);
        end

        function functional = createContinuumDamageFunctional(obj)
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.material           = obj.createDamagedMaterial();
            s.quadOrder          = 2;
            s.test               = LagrangianFunction.create(obj.mesh,2,'P1');
            functional = ContinuumDamageFunctional(s);
        end        
    
        function dM = createDamagedMaterial(obj)
            s.type = 'ContinuumDamage';
            s.baseMaterial = obj.createBaseMaterial();
            s.damage       = obj.createDamagedLaw();
            dM = Material.create(s);
        end

        function mat = createBaseMaterial(obj)
            E  = obj.matInfo.young;
            nu = obj.matInfo.poisson;
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
            r0     = obj.damageInfo.r0;
            s.r0   = ConstantFunction.create(r0,obj.mesh);
            s.type = obj.damageInfo.type;   
            s.params = obj.damageInfo.params;
            hL = HardeningLaw.create(s);
        end          
 
    end
end