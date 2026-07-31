classdef Tutorial08_ContinuumDamage < handle

    properties (Access = public)
        output
    end

    properties (Access = private)
        mesh
        boundaryConditions
        r0, internalDamageVariable
        functional
        damage
    end

    methods (Access = public)

        function obj = Tutorial08_ContinuumDamage()
            obj.init();
            obj.defineCase();
            obj.createInternalDamageVariable();
            obj.createContinuumDamageFunctional();
            obj.solveContinuumDamageProblem();
        end

        function solveContinuumDamageProblem(obj)
            s.mesh                   = obj.mesh;
            s.boundaryConditions     = obj.boundaryConditions;
            s.internalDamageVariable = obj.internalDamageVariable;
            s.functional             = obj.functional;

            s.monitoring.set         = true;
            s.monitoring.print       = true;

            s.tolerance              = 1e-8;
            s.maxIter                = 20;

            CDComp = ContinuumDamageComputer(s);
            obj.output = CDComp.compute();
        end
    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function defineCase(obj)
            s.mesh.type   = 'Rectangle';
            s.mesh.length = 1;
            s.mesh.width  = 1;
            s.mesh.lN     = 1;
            s.mesh.wN     = 1;
            s.bc.type   = 'DisplacementTractionY';
            s.bc.values = [0:1e-1:2];
            [obj.mesh, obj.boundaryConditions] = BenchmarkManager.create(s);
        end

        function createInternalDamageVariable(obj)
            obj.r0 = 10;
            s.r0   = ConstantFunction.create(obj.r0,obj.mesh);
            s.mesh = obj.mesh;
            obj.internalDamageVariable = InternalDamageVariable(s);
        end

        function createContinuumDamageFunctional(obj)
            [bMat, sMat,tMat]    = obj.createDamagedMaterials();
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.bMat               = bMat;
            s.sFun               = sMat;
            s.tFun               = tMat;
            s.quadOrder          = 2;
            s.test               = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.functional       = ContinuumDamageFunctional(s);
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
            s.r0   = ConstantFunction.create(obj.r0,obj.mesh);
            s.type = 'Linear';
            s.params.r1 = 20;
            s.params.hardening = -0.5;
            % s.params.A = 1;
            % s.params.qInf = 0;
            % s.params.w = 500;
            hL = HardeningLaw.create(s);
        end          
 
    end
end