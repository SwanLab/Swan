classdef TestingContinuumDamage < handle

    properties (Access = public)
        data
    end

    properties (Access = private)
        tol
    end

    properties (Access = private)
        mesh
        bc
        internalDamageVariable
        damageFunctional
    end

    methods (Access = public)

        function obj = TestingContinuumDamage()
            obj.init();
            obj.mesh             = obj.createMesh();
            obj.bc               = obj.createBoundaryConditions();
            obj.damageFunctional = obj.createContinuumDamageFunctional();
            obj.internalDamageVariable = obj.createInternalDamageVariable();
            c = obj.createContinumDamageComputer();
            obj.data = c.compute();
        end
    end

    methods (Access = private)

        function init(obj)
            close all
            obj.tol = 1e-6;
        end

        function mesh = createMesh(obj)
                l = 1;
                w = 1;
                N = 40;
                M = 40;
                mesh = QuadMesh(l,w,N,M);

                % file = 'SENshear0_0025';
                % a.fileName = file;
                % s = FemDataContainer(a);
                % mesh = s.mesh;
        end

        function bc = createBoundaryConditions(obj)
            s.mesh = obj.mesh;
            s.bcType = 'displacementMixed'; %fiberMatrix
            s.bcValueSet = [0:1e-3:2];            
            bc = BcContinuumDamage(s);
        end

        function r = createInternalDamageVariable(obj)
            s.r0 = obj.createR0();
            s.mesh = obj.mesh;
            r = InternalDamageVariable(s);
        end

        function sF = createContinuumDamageFunctional(obj)
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.bc;
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
            r1     = 25;
            s.r1   = ConstantFunction.create(r1,obj.mesh);
            s.type = 'Exp';%'Exp'; %'Exp'            
            s.H    = -0.1;
            s.A    = 0.2;
            s.r0   = obj.createR0();
            s.w1   = 80;
            s.qInf = -1;
            hL = HardeningLaw.create(s);
        end          

        function r0 = createR0(obj)
            r0 = 10;
            r0 = ConstantFunction.create(r0,obj.mesh);            
        end

        function comp = createContinumDamageComputer(obj)
            s.mesh               = obj.mesh;
            s.boundaryConditions = obj.bc;
            s.tol                = obj.tol;
            s.damageFunctional   = obj.damageFunctional();
            s.internalDamageVariable = obj.internalDamageVariable;
            comp = ContinuumDamageComputer(s);
        end        
 
    end
end