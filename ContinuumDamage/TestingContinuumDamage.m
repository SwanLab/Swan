classdef TestingContinuumDamage < handle
    
    properties (Access = public)
        mesh
        bc
        material
        results
    end
    
    methods (Access = public)
        
        function obj = TestingContinuumDamage(cParams)
            obj.mesh     = obj.createMesh(cParams);
            obj.bc       = obj.createBoundaryConditions(cParams);
            obj.material = obj.createMaterial(cParams);
            obj.results  = obj.compute(cParams);
        end

    function compareWithElasticProblem(obj,dataIn)
            s.mesh = obj.mesh;
            s.material = obj.material;
            s.boundaryConditions = obj.bc;
            s.type = dataIn.type;
            s.scale = dataIn.scale;
            s.solverType = dataIn.solverType;
            s.solverMode = dataIn.solverMode;
            s.solverCase = dataIn.solverCase;

            EP = ElasticProblem(s);
            Ref = EP.solve();
           
            if (Ref == obj.results)

                fprintf ("Continuum Damage TEST: \nPASSED\n") %Optional message
                disp ("-------------------")
            else
                disp ("Continuum Damage TEST:")
                fprintf (2,'FAILED\n')
                disp ("-------------------")
            end
    end
end
    
    methods (Access = private)
        
        function mesh = createMesh(obj,s)
            l = s.meshLength;
            w = s.meshWidth;
            N = s.meshN;
            M = s.meshM;
            mesh = QuadMesh(l,w,N,M);
        end
        
        function bc = createBoundaryConditions(obj,s)
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
        end
        
        function mat = createMaterial(obj,s)
            sMat.ndim    = s.ndim;
            sMat.young   = AnalyticalFunction.create(@(x) s.E*ones(size(x,[2,3])),1,obj.mesh);
            sMat.poisson = AnalyticalFunction.create(@(x) s.nu*ones(size(x,[2,3])),1,obj.mesh);
            mat = Isotropic2dElasticMaterial(sMat);
        end
        
        function results = compute(obj,sSolver)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.bc;
            s.material = obj.material;
            s.type = sSolver.type;
            s.scale = sSolver.scale;
            s.solverType = sSolver.solverType;
            s.solverMode = sSolver.solverMode;
            s.solverCase = sSolver.solverCase;

            comp = ContinuumDamageComputer(s);
            results = comp.compute();
        end
    end
end