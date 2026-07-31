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
        damage
    end

    methods (Access = public)

        function obj = TestingContinuumDamage(cParams)
            obj.init(cParams);
            [obj.mesh, obj.boundaryConditions] = obj.defineCase();
            obj.internalDamageVariable = obj.createInternalDamageVariable();
            obj.functional             = obj.createContinuumDamageFunctional();
            
        end

        function outputData = compute(obj)
            s.mesh                   = obj.mesh;
            s.boundaryConditions     = obj.boundaryConditions;
            s.internalDamageVariable = obj.internalDamageVariable;
            s.functional             = obj.functional;
            s.damage                 = obj.damage;
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

        function [mesh,bc] = defineCase(obj)
            [mesh, bc] = BenchmarkManager.create(obj.benchmark);
        end

        function r = createInternalDamageVariable(obj)
            r0     = obj.damageInfo.r0;
            s.r0   = ConstantFunction.create(r0,obj.mesh);
            s.mesh = obj.mesh;
            r = InternalDamageVariable(s);
        end

        function functional = createContinuumDamageFunctional(obj)
            [bMat, sMat,tMat]    = obj.createDamagedMaterials();
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.bMat               = bMat;
            s.sMat               = sMat;
            s.tMat               = tMat;
            s.quadOrder          = 2;
            s.test               = LagrangianFunction.create(obj.mesh,2,'P1');
            functional       = ContinuumDamageFunctional(s);
     
        end        
    
        function [baseMat, secantMat, tangentMat] = createDamagedMaterials(obj)
            baseMat    = obj.createBaseMaterial();
            d          = obj.createDamagedLaw();
            secantMat  = @(r) ContinuumDamageMaterials.obtainTensorSecant(baseMat,d,r);
            tangentMat = @(u,r) ContinuumDamageMaterials.obtainTensorTangent(baseMat,d,u,r);
        end

        function mat = createBaseMaterial(obj)
            N      = obj.mesh.ndim;
            E      = ConstantFunction.create(210,obj.mesh);
            nu     = ConstantFunction.create(0.3,obj.mesh);
            mu     = E./(2*(1+nu));  mu = Expand(mu,4);
            lambda = LameLambda(E,nu,N);  lambda = Expand(lambda,4);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            mat    = 2*mu.*I + lambda.*IxI;
        end


        function d = createDamagedLaw(obj)
            s.hardeningLaw = obj.createHardeningLaw();
            d = DamageLaw(s);
            obj.damage = d;
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