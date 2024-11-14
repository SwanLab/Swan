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
            l = s.meshLength;
            w = s.meshWidth;
            N = s.meshN;
            M = s.meshM;
            mesh = QuadMesh(l,w,N,M);
        end

        function bc = defineBoundaryConditions(obj,s)
            switch s.bcType
                case 'displacementTraction'
                    isInDown = @(coord) (abs(coord(:,2) - min(coord(:,2)))< 1e-12);
                    sDir.domain    = @(coor) isInDown(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isInUp = @(coord) (abs(coord(:,2) - max(coord(:,2)))< 1e-12);
                    sDir.domain    = @(coor) isInUp(coor);
                    sDir.direction = [2];
                    sDir.value     = s.bcVal;
                    Dir2 = DirichletCondition(obj.mesh,sDir);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1 Dir2];
                    s.pointloadFun = [];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);

                case 'forceTraction'
                    isInDown = @(coord) (abs(coord(:,2) - min(coord(:,2)))< 1e-12);
                    sDir.domain    = @(coor) isInDown(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isInUp = @(coord) (abs(coord(:,2) - max(coord(:,2)))< 1e-12);
                    sNeum.domain    = @(coor) isInUp(coor);
                    sNeum.direction = [2];
                    sNeum.value     = s.bcVal;
                    Neum1 = PointLoad(obj.mesh,sNeum);


                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1];
                    s.pointloadFun = [Neum1];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);

                case 'displacementBending'
                    isInLeft = @(coord) (abs(coord(:,1) - min(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isInLeft(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isInRight = @(coord) (abs(coord(:,1) - max(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isInRight(coor);
                    sDir.direction = [2];
                    sDir.value     = s.bcVal;
                    Dir2 = DirichletCondition(obj.mesh,sDir);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1 Dir2];
                    s.pointloadFun = [];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);

                case 'forceBending'
                    isInLeft = @(coord) (abs(coord(:,1) - min(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isInLeft(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                    Dir1 = DirichletCondition(obj.mesh,sDir);

                    isInRight = @(coord) (abs(coord(:,1) - max(coord(:,1)))< 1e-12);
                    sDir.domain    = @(coor) isInRight(coor);
                    sDir.direction = [2];
                    sDir.value     = s.bcVal;
                    Dir2 = PointLoad(obj.mesh,sDir);

                    s.mesh = obj.mesh;
                    s.dirichletFun = [Dir1 Dir2];
                    s.pointloadFun = [];
                    s.periodicFun = [];
                    bc = BoundaryConditions(s);  
            end
        end

        function mat = createMaterial(obj,s)
         
            s.mesh = obj.mesh;
            
            mat = DamagedMaterial(s);
        end

    end
end