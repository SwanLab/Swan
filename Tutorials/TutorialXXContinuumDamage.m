classdef TutorialXXContinuumDamage < handle

    properties (Access = public)
        output
    end

    properties (Access = private)
        mesh
        boundaryConditions
        r0, internalDamageVariable
        functional
    end

    methods (Access = public)

        function obj = TutorialXXContinuumDamage()
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
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.material           = obj.createDamagedMaterial();
            s.quadOrder          = 2;
            s.test               = LagrangianFunction.create(obj.mesh,2,'P1');
            obj.functional = ContinuumDamageFunctional(s);
        end        
    
        function dM = createDamagedMaterial(obj)
            s.type = 'ContinuumDamage';
            s.baseMaterial = obj.createBaseMaterial();
            s.damage       = obj.createDamagedLaw();
            dM = Material.create(s);
        end

        function mat = createBaseMaterial(obj)
            E  = 210;
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