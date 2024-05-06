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

        neohookeanFun

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
            close all;
            obj.init();
            obj.createDisplacementFun();
            obj.createBoundaryConditions();
            obj.computeForces();

            s.material = obj.material;
            s.mesh = obj.mesh;
            neo = NeohookeanFunctional(s);
            obj.neohookeanFun = neo;
%             fint = neo.computeInternalForces(obj.uFun);
%             hess = neo.computeHessian(obj.uFun);
            
            

            % Check first piola convergence

            bc = obj.boundaryConditions;
            u_k = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;

            r = 1;
            i = 1;
            rpre = 1;
            alpha = 0.01;
            f = animatedline;
            obj.applyDirichletToUFun();
            while r > 10e-6
                val = max(neo.compute(obj.uFun))
                Fint = obj.computeInternalForces();
                hess = neo.computeHessian(obj.uFun);
                res  = Fint - obj.Fext;
%                 deltaUk = hess\res;
                u_next = u_k - alpha*res;
%                 u_next = u_k - deltaUk;
                u_next(bc.dirichlet_dofs) = bc.dirichlet_vals;
                obj.uFun.fValues = reshape(u_next,[obj.mesh.ndim,obj.mesh.nnodes])';
                r = norm(u_next - u_k)
                u_k = u_next;
                i = i+1;
%                 if r>1
%                     break
%                 end
                addpoints(f,i,r);
                drawnow
                rpre = r;
            end
% 
% 
%             bc = obj.boundaryConditions;
%             xpre = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
%             xpre(bc.dirichlet_dofs) = bc.dirichlet_vals;
%             nIter = 0;
%             while nIter < 10
%                 Fint = obj.computeInternalForces();
%                 res  = Fint + obj.Fext;
%                 hess = obj.computeSecondPiola();
%                 [hess_red,res_red] = obj.full2reduced(hess,res);
%                 x = hess_red\res_red;
%                 % x = obj.computeNewtonRaphson(xpre, res, hess);
%                 obj.uFun.fValues = reshape(x,[obj.mesh.ndim,obj.mesh.nnodes])';
%                 norm(x-xpre)
%                 xpre = x;
%                 nIter = nIter + 1;
%             end
        end

        function solve(obj)
        end

    end

    methods (Access = private)

        function init(obj)
%             obj.mesh = UnitHexaMesh(5,5,5);
            obj.mesh = UnitQuadMesh(1,1);
            obj.material.lambda = 3/4;
            obj.material.mu = 3/8;
        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            obj.uFun.fValues = obj.uFun.fValues + 0;
        end

        function applyDirichletToUFun(obj)
            bc = obj.boundaryConditions;
            u_k = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;
            obj.uFun.fValues = reshape(u_k,[obj.mesh.ndim,obj.mesh.nnodes])';
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
%             s.mesh = obj.mesh;
%             s.material = obj.material;
%             test = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
%             rhs = RHSintegrator_FirstPiola(s);
%             intfor = rhs.compute(obj.uFun, test);
            intfor = obj.neohookeanFun.computeInternalForces(obj.uFun);
        end

        function hess = computeSecondPiola(obj)
            s.mesh = obj.mesh;
            s.test  = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            s.trial = obj.uFun;
            s.material = obj.material;
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
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight = @(coor)  abs(coor(:,1))==xMax;

            sDir.domain    = @(coor) isLeft(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            s.dirichletFun =  DirichletCondition(obj.mesh, sDir);


            sPL.domain    = @(coor) isRight(coor);
            sPL.direction = 2;
            sPL.value     = -1;
            s.pointloadFun = PointLoad(obj.mesh, sPL);


%             sDir{1}.domain    = @(coor) isLeft(coor);
%             sDir{1}.direction = 1;
%             sDir{1}.value     = -0.1;
%             s.dirichletFun =  DirichletCondition(obj.mesh, sDir);
% 
% 
%             sDir{2}.domain    = @(coor) isLeft(coor);
%             sDir{2}.direction = [2,3];
%             sDir{2}.value     = 0;
%             s.dirichletFun =  DirichletCondition(obj.mesh, sDir);
% 
% 
%             sDir{3}.domain    = @(coor) isRight(coor);
%             sDir{3}.direction = 1;
%             sDir{3}.value     = +0.1;
%             s.dirichletFun =  DirichletCondition(obj.mesh, sDir);
% 
% 
%             sDir{4}.domain    = @(coor) isRight(coor);
%             sDir{4}.direction = [2,3];
%             sDir{4}.value     = 0;
%             s.dirichletFun =  DirichletCondition(obj.mesh, sDir);

            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

        function [LHS,RHS] = full2reduced(obj,LHS,RHS)
            bc = obj.boundaryConditions;
            dofs = 1:obj.uFun.nDofs;
            free_dofs = setdiff(dofs, bc.dirichlet_dofs);
            LHS = LHS(free_dofs, free_dofs);
            RHS = RHS(free_dofs);
        end

    end

end
