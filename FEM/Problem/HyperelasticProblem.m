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
%             fint = obj.computeIntForcesShape(perc);
            hess = neo.computeHessian(obj.uFun);
            
            

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

            nsteps = 5;
            F = zeros(obj.uFun.nDofs,1);
            R = zeros(obj.uFun.nDofs,1);
            for iStep = 1:nsteps
                disp('------')
                disp('------')
                disp('------')
                disp('NEW LOAD STEP')
                loadPercent = iStep/nsteps;
%                 DeltaF = -obj.computeForces(loadPercent);
                DeltaF = obj.computeExtForcesShape(loadPercent);
                F = F + DeltaF;
                R = R - DeltaF;
                while norm(R)/norm(F) > 10e-8
                    % Energy
                    val = max(neo.compute(obj.uFun))
    
                    % Hessian
                    hess = neo.computeHessian(obj.uFun);
                    h_red = hess(free,free);

                    deltaUk_free = -h_red\R(free);
                    deltaUk = zeros(size(R));
                    deltaUk(free) = deltaUk_free;
%                     u_next = u_k + deltaUk;
                    u_next = deltaUk;
    
%                     u_next = u_k - alpha*res;
                    u_next(bc.dirichlet_dofs) = bc.dirichlet_vals;
                    obj.uFun.fValues = reshape(u_next,[obj.mesh.ndim,obj.mesh.nnodes])';


                    sigma = obj.computeCauchyStress();
    
                    % Residual
                    T = obj.computeInternalForces();
                    R  = T - F;

                    r = norm(R)/norm(F)
                    u_k = u_next;
                    i = i+1;
                    addpoints(f,i,r);
                    drawnow
                end
            end
% 
        end

        function solve(obj)
        end

    end

    methods (Access = private)

        function init(obj)
%             obj.mesh = HexaMesh(2,1,1,20,5,5);
%             obj.mesh = UnitHexaMesh(5,5,5);
            obj.mesh = UnitQuadMesh(6,6);
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

        function Fext = computeExtForcesShape(obj,perc)
            pl = obj.boundaryConditions.pointloadFun;
            pl.fValues = pl.fValues*perc;
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = 1;
            rhsI       = RHSintegrator.create(s);
            test = LagrangianFunction.create(obj.mesh,obj.uFun.ndimf,'P1');
            Fext = rhsI.compute(pl,test);
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

        function sigma = computeCauchyStress(obj)
            mu = obj.material.mu;
            lambda = obj.material.lambda;
            GradU2 = ActualGrad(obj.uFun);
            Id = Identity(obj.uFun);
            F = GradU2 + Id;
            b = F*F';
            jac = Det(F);
            sigma = mu.*(b-Id)./jac + lambda.*(log(jac)).*Id./jac;
        end

        function [F,I33] = computeDeformationGradient(obj, uFun, xG)
            nPoints  = size(xG,2);
            nElem = obj.mesh.nelem;
            nDimG = obj.mesh.ndim;
            nDimf = uFun.ndimf;

            GradU = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);
            GradU = permute(GradU, [2 1 3 4]);
%             GradU = [0.0 0.0 0.0; -3.415063509461096 -0.24999999999999956 -0.4330127018922192; 0.9150635094610968 0.43301270189221924 -0.24999999999999994];

            I33 = obj.createIdentityMatrix(size(GradU));

            F = I33 + GradU;
%             if size(F,1) == 2
%                 F = [F,zeros(2,1,nPoints,nElem); zeros(1,2,nPoints,nElem), ones(1,1,nPoints,nElem)];
%                 I33 = obj.createIdentityMatrix(size(F));
%             end
        end
        
        function createBoundaryConditions(obj)
%             obj.createBC2D_oneelem();
            obj.createBC2D_nelem();
%             obj.createBC3D_oneelem();
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
            sPL.value     = 0.1;
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
            
            % 3D N ELEMENTS
%             sDir1.domain    = @(coor) isLeft(coor) & ~isMiddle(coor);
%             sDir1.direction = [1];
%             sDir1.value     = 0;
%             dir1 =  DirichletCondition(obj.mesh, sDir1);
% 
%             sDir2.domain    = @(coor) isLeft(coor) & isMiddle(coor);
%             sDir2.direction = [1,2,3];
%             sDir2.value     = 0;
%             dir2 =  DirichletCondition(obj.mesh, sDir2);
%             s.dirichletFun = [dir1, dir2];

            sDir2.domain    = @(coor) isLeft(coor);
            sDir2.direction = [1,2,3];
            sDir2.value     = 0;
            s.dirichletFun =  DirichletCondition(obj.mesh, sDir2);

            sPL.domain    = @(coor) isRight(coor);
            sPL.direction = 1;
            sPL.value     = 0.0000001;
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
