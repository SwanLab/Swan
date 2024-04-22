classdef HyperelasticProblem < handle
    
    properties (Access = public)
        uFun
        strainFun
        stressFun
    end

    properties (Access = private)
        quadrature
        boundaryConditions, BCApplier
        dirichletFun

        stiffness
        Fext
        solver, solverType, solverMode, solverCase
        scale
        
        strain, stress
    end

    properties (Access = protected)
        mesh 
        material  
        displacementFun
    end

    methods (Access = public)

        function obj = HyperelasticProblem()
            obj.init();
            obj.createDisplacementFun();
            obj.createBoundaryConditions();
            obj.computeForces();
            bc = obj.boundaryConditions;
            xpre = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
            xpre(bc.dirichlet_dofs) = bc.dirichlet_vals;
            nIter = 0;
            while nIter < 10
                Fint = obj.computeInternalForces();
                res  = Fint - obj.Fext;
                hess = obj.computeSecondPiola();
                x = obj.computeNewtonRaphson(xpre, res, hess);
                obj.uFun.fValues = reshape(x,[obj.mesh.ndim,obj.mesh.nnodes])';
                norm(x-xpre)
                xpre = x;
                nIter = nIter + 1;
            end
        end

        function solve(obj)
        end

    end

    methods (Access = private)

        function init(obj)
            obj.mesh = UnitHexaMesh(3,3,3);
            obj.material.lambda = 1;
            obj.material.mu = 1;
        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            obj.uFun.fValues = obj.uFun.fValues + 0.001;
        end

        function computeForces(obj)
            s.type     = 'Elastic';
            s.scale    = 'MACRO';
            s.dim.ndofs = obj.uFun.nDofs;
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            obj.Fext = rhs;
        end

        function intfor = computeInternalForces(obj)
            s.mesh = obj.mesh;
            test = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            rhs = RHSintegrator_FirstPiola(s);
            intfor = rhs.compute(obj.uFun, test);
        end

        function hess = computeSecondPiola(obj)
            s.mesh = obj.mesh;
            s.test  = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            s.trial = obj.uFun;
            lhs = LHSintegrator_SecondPiola(s);
            hess = lhs.compute();
        end

        function x = computeNewtonRaphson(obj, xpre, res, hess)
            bc = obj.boundaryConditions;
            h = 1./hess;
            h = diag(h);
            h(isinf(h)) = 1;
            x = xpre - h.*res;
            x(bc.dirichlet_dofs) = bc.dirichlet_vals;
        end


        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  abs(coor(:,1))==xMax;

            sDir.domain    = @(coor) isDir(coor);
            sDir.direction = [1,2,3];
            sDir.value     = 0;
            s.dirichletFun =  DirichletCondition(obj.mesh, sDir);


            sPL.domain    = @(coor) isForce(coor);
            sPL.direction = 2;
            sPL.value     = -1;
            s.pointloadFun = PointLoad(obj.mesh, sPL);

            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end
        
    end

end
