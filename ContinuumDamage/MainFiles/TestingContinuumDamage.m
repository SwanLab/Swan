classdef TestingContinuumDamage < handle
    properties (Access = private)
        mesh
        bc
        material
        solverParams
        qLaw

        H
        r0
        r1
        tol
    end

    methods (Access = public)
        function obj = TestingContinuumDamage(cParams)
            obj.mesh         = obj.createMesh(cParams.mesh);
            obj.bc           = obj.defineBoundaryConditions(cParams.bc);
            obj.material     = obj.createMaterial(cParams.material);
            obj.solverParams = cParams.solver;
            obj.H            = cParams.H;
            obj.r0           = cParams.r0;
            obj.r1           = cParams.r1;
            obj.tol          = cParams.tol;
            obj.qLaw         = cParams.qLaw;
        end

        function data = compute(obj)
            sComp.mesh = obj.mesh;
            sComp.boundaryConditions = obj.bc;
            sComp.material = obj.material;
            sComp.solver = obj.solverParams;
            sComp.H = obj.H;
            sComp.r0 = obj.r0;
            sComp.r1 = obj.r1;
            sComp.tol = obj.tol;
            sComp.qLaw = obj.qLaw;
            comp = ContinuumDamageComputer(sComp);
            data = comp.compute();
        end

        function compareWithElasticProblem(~,data,uRef)
            if  all(all(ismembertol(uRef.fValues,data.displacement.fValues,1e-10)))
                fprintf ("Continuum Damage TEST: \nPASSED\n")
                disp ("-------------------")
            else
                disp ("Continuum Damage TEST:")
                fprintf (2,'FAILED\n')
                disp ("-------------------")
            end
        end
    end

    methods (Access = private)
        function mesh = createMesh(~,s)
            if ~isfield(s,'name')
                l = s.meshLength;
                w = s.meshWidth;
                N = s.meshN;
                M = s.meshM;
                mesh = QuadMesh(l,w,N,M);
            else
                file = s.name;
                a.fileName = file;
                s = FemDataContainer(a);
                mesh = s.mesh;
            end
        end

        function bc = defineBoundaryConditions(obj,s)
            s.mesh = obj.mesh;
            bc = BcContinuumDamage(s);
        end

        function mat = createMaterial(obj,s)
            s.mesh = obj.mesh;
            mat = DamagedMaterial(s);
        end
    end
end