classdef TestingContinuumDamage < handle

    properties (Access = private)
        mesh
        bc
        material
        solverParams

    end

    methods (Access = public)

        function obj = TestingContinuumDamage(cParams)
            obj.mesh      = obj.createMesh(cParams.mesh);
            obj.bc        = obj.defineBoundaryConditions(cParams.bc);
            obj.material  = obj.createMaterial(cParams.material);
            obj.solverParams = cParams.solver;


        end

        function data = compute(obj)
            sComp.mesh = obj.mesh;
            sComp.boundaryConditions = obj.bc;
            sComp.material = obj.material;
            sComp.solver = obj.solverParams;

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
            bc = Bc_ContinuumDamage(s);
        end

        function mat = createMaterial(obj,s)
         
            s.mesh = obj.mesh;
            
            mat = DamagedMaterial(s);
        end

    end
end