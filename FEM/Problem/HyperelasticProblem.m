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

            s.material = obj.material;
            s.mesh = obj.mesh;
            neo = NeohookeanFunctional(s);
            obj.neohookeanFun = neo;
%             fint = neo.computeInternalForces(obj.uFun);
%             hess = neo.computeHessian(obj.uFun);
            
            

            % Check first piola convergence

            bc = obj.boundaryConditions;
            dofs = 1:obj.uFun.nDofs;
            free = setdiff(dofs, bc.dirichlet_dofs);
            u_k = reshape(obj.uFun.fValues',[obj.uFun.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;

            r = 1;
            i = 1;
            rpre = 1;
            alpha = 0.001;
            f = animatedline;
            obj.applyDirichletToUFun();

            nsteps = 1;
            for iStep = 1:nsteps
                disp('------')
                disp('------')
                disp('------')
                disp('NEW LOAD STEP')
                loadPercent = iStep/nsteps;
                Fext = obj.computeForces(loadPercent);
                while r > 10e-8
                    % Energy
                    val = max(neo.compute(obj.uFun))
    
                    % Residual
                    Fint = obj.computeInternalForces();
                    res  = Fint - Fext;
    
%                     % Hessian
%                     hess = neo.computeHessian(obj.uFun);
%                     h_red = hess(free,free);
% 
%                     deltaUk_free = h_red\res(free);
%                     deltaUk = zeros(size(Fint));
%                     deltaUk(free) = deltaUk_free;
%                     u_next = u_k - deltaUk;
    
                    u_next = u_k - alpha*res;
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
%             obj.mesh = HexaMesh(2,1,1,20,5,5);
            obj.mesh = UnitHexaMesh(1,1,1);
%             obj.mesh = UnitQuadMesh(5,5);
%             obj.material.lambda = 3/4;
%             obj.material.mu = 3/8;
            E = 10.0;
            nu = 0.3;
            obj.material.lambda = E*nu/((1 + nu)*(1 - 2*nu));
            obj.material.mu = E/(2*(1 + nu));
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

        function Fext = computeForces(obj, perc)
            s.type     = 'Elastic';
            s.scale    = 'MACRO';
            s.dim.ndofs = obj.uFun.nDofs;
            s.BC       = obj.boundaryConditions;
            s.BC.pointload_vals = s.BC.pointload_vals*perc;
            s.mesh     = obj.mesh;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            Fext = rhs;
        end

        function intfor = computeInternalForces(obj)
            intfor = obj.neohookeanFun.computeInternalForces(obj.uFun,obj.boundaryConditions);
        end

        function createBoundaryConditions(obj)
%             obj.createBC2D_oneelem();
%             obj.createBC2D_nelem();
            obj.createBC3D_oneelem();
%             obj.createBC3D_nelem();
        end

        function bc = createBC2D_oneelem(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isTop    = @(coor)  abs(coor(:,2))==yMax;
            isBottom = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(coor(:,2))==yMax/2;

            % 2D ONE ELEMENT
            sDir1.domain    = @(coor) isLeft(coor) & isTop(coor);
            sDir1.direction = [1];
            sDir1.value     = 0;
            dir1 =  DirichletCondition(obj.mesh, sDir1);

            sDir2.domain    = @(coor) isLeft(coor) & isBottom(coor);
            sDir2.direction = [1,2];
            sDir2.value     = 0;
            dir2 =  DirichletCondition(obj.mesh, sDir2);
            s.dirichletFun = [dir1, dir2];

            sPL.domain    = @(coor) isRight(coor);
            sPL.direction = 1;
            sPL.value     = 1;
            s.pointloadFun = PointLoad(obj.mesh, sPL);
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

        function bc = createBC2D_nelem(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isTop    = @(coor)  abs(coor(:,2))==yMax;
            isBottom = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(coor(:,2))==yMax/2;
            
            % 2D N ELEMENTS
            sDir1.domain    = @(coor) isLeft(coor) & ~isMiddle(coor);
            sDir1.direction = [1];
            sDir1.value     = 0;
            dir1 =  DirichletCondition(obj.mesh, sDir1);

            sDir2.domain    = @(coor) isLeft(coor) & isMiddle(coor);
            sDir2.direction = [1,2];
            sDir2.value     = 0;
            dir2 =  DirichletCondition(obj.mesh, sDir2);
            s.dirichletFun = [dir1, dir2];

            sPL.domain    = @(coor) isRight(coor);
            sPL.direction = 1;
            sPL.value     = 1;
            s.pointloadFun = PointLoad(obj.mesh, sPL);
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

        function bc = createBC3D_oneelem(obj)
            xMax    = max(obj.mesh.coord(:,1));
            zMax    = max(obj.mesh.coord(:,3));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isTop    = @(coor)  abs(coor(:,3))==zMax;
            isBottom = @(coor)  abs(coor(:,3))==0;
            isMiddle = @(coor)  abs(coor(:,3))==zMax/2;

            % 2D ONE ELEMENT
            sDir1.domain    = @(coor) isLeft(coor) & isTop(coor);
            sDir1.direction = [1];
            sDir1.value     = 0;
            dir1 =  DirichletCondition(obj.mesh, sDir1);

            sDir2.domain    = @(coor) isLeft(coor) & isBottom(coor);
            sDir2.direction = [1,2,3];
            sDir2.value     = 0;
            dir2 =  DirichletCondition(obj.mesh, sDir2);
            s.dirichletFun = [dir1, dir2];

            sPL.domain    = @(coor) isRight(coor);
            sPL.direction = 1;
            sPL.value     = 1;
            s.pointloadFun = PointLoad(obj.mesh, sPL);
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

        function bc = createBC3D_nelem(obj)
            xMax    = max(obj.mesh.coord(:,1));
            zMax    = max(obj.mesh.coord(:,3));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isTop    = @(coor)  abs(coor(:,3))==zMax;
            isBottom = @(coor)  abs(coor(:,3))==0;
            isMiddle = @(coor)  abs(coor(:,3))==zMax/2;
            
            % 2D N ELEMENTS
            sDir1.domain    = @(coor) isLeft(coor) & ~isMiddle(coor);
            sDir1.direction = [1];
            sDir1.value     = 0;
            dir1 =  DirichletCondition(obj.mesh, sDir1);

            sDir2.domain    = @(coor) isLeft(coor) & isMiddle(coor);
            sDir2.direction = [1,2,3];
            sDir2.value     = 0;
            dir2 =  DirichletCondition(obj.mesh, sDir2);
            s.dirichletFun = [dir1, dir2];

            sPL.domain    = @(coor) isRight(coor);
            sPL.direction = 1;
            sPL.value     = 1;
            s.pointloadFun = PointLoad(obj.mesh, sPL);
            
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
