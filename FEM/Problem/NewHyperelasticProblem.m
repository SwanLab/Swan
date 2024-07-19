classdef NewHyperelasticProblem < handle
    
    properties (Access = public)
        uFun
        strainFun
        stressFun
    end

    properties (Access = private)
        quadrature
        boundaryConditions, BCApplier
        dirichletFun

        neohookeanFun, linearFun

        stiffness
        %Fext
        solver, solverType, solverMode, solverCase
        scale
        
        strain, stress

        FextInitial
    end

    properties (Access = protected)
        mesh 
        material, materialLinear
        displacementFun
    end

    methods (Access = public)

        function obj = NewHyperelasticProblem()
            obj.init();
            u = obj.uFun;
            
            nSteps = 10;

            for iStep = 1:nSteps
                % Calculate boundary conditions
                perc = iStep/nSteps;
                obj.createBoundaryConditions(perc);
                obj.applyDirichletToUFun(u);

                hasNotConverged = 1;
                while hasNotConverged
                    LHS = obj.computeLHS(u);
                    RHS = obj.computeResidual(u);
                    u = obj.computeDisplacements();
%                     F = obj.computeForces();
%                     Fint = obj.computeInternalForces();
                end
            end

        end

    end

    methods (Access = private)

        function init(obj)
            obj.createMesh();
            obj.createMaterial();
            obj.createDisplacement();
            obj.createFunctionals();
        end

        function createMesh(obj)
%             obj.mesh = HexaMesh(2,1,1,20,5,5);
%             obj.mesh = UnitHexaMesh(15,15,15);
%             obj.mesh = UnitQuadMesh(14,14);
            IM = Mesh.createFromGiD('hole_mesh_quad.m');
            obj.mesh = IM;
        end
        
        function createMaterial(obj)
            % Only 3D
            obj.material.mu = 1*1000;        % kPa
            obj.material.lambda = 1*10*1000; % kPa
            obj.createLinearMaterial();
        end


        function createLinearMaterial(obj)
            G = obj.material.mu;
            L = obj.material.lambda;

            E1 = G*(3*L+2*G)/(L+G);
            nu1 = L / (2*(L+G));
            E         = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            nu        = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,obj.mesh);
            s.pdim    = obj.mesh.ndim;
            s.nelem   = obj.mesh.nelem;
            s.mesh    = obj.mesh;
            s.young   = E;
            s.poisson = nu;
            s.type = 'ISOTROPIC';
            s.ndim = obj.mesh.ndim;
            mat = Material.create(s);
            obj.materialLinear = mat;
        end

        function createDisplacement(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function createFunctionals(obj)
            s.material = obj.material;
            s.mesh     = obj.mesh;
            neo = NeohookeanFunctional(s);
            obj.neohookeanFun = neo;

            s.lambda = s.material.lambda;
            s.mu = s.material.mu;
            s.material = obj.materialLinear;
            lin = LinearElasticityFunctional(s);
            obj.linearFun = lin;
        end


        %% Loop methods

        function intfor = computeInternalForces(obj)
            intfor = obj.neohookeanFun.computeInternalForces(obj.uFun,obj.boundaryConditions);
        end

        function LHS = computeLHS(obj, uFun)
            LHS = obj.neohookeanFun.computeHessian(uFun);
            LHS = obj.linearFun.computeHessian(uFun);
        end

        function RHS = computeResidual(obj,uFun)
            Fint = obj.neohookeanFun.computeGradient(uFun);

            Fint2 = obj.linearFun.computeGradient(uFun);
            RHS = Fint;
        end

        function uFun = computeDisplacements(obj)
        end


        function [Fext] = computeExternalForces(obj,perc)
            Fext = perc*obj.FextInitial;
            Fext = obj.reshapeToVector(Fext);
        end

        %% Boundary conditions

        function applyDirichletToUFun(obj, uFun)
            bc = obj.boundaryConditions;
            u_k = reshape(uFun.fValues',[uFun.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;
            uFun.fValues = reshape(u_k,[obj.mesh.ndim,obj.mesh.nnodes])';
        end
        
        function createBoundaryConditions(obj,perc)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft   = @(coor)  abs(coor(:,1))==0;
            isRight  = @(coor)  abs(coor(:,1))==xMax;
            isHalf   = @(coor)  abs(abs(coor(:,1)) - xMax/2) <=10e-2;
            isTop    = @(coor)  abs(coor(:,2))==yMax;
            isBottom = @(coor)  abs(coor(:,2))==0;
            isMiddle = @(coor)  abs(abs(coor(:,2))-yMax/2) <= 10e-2;
            
            % 2D N ELEMENTS
            sDir.domain    = @(coor) isBottom(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            dir =  DirichletCondition(obj.mesh, sDir);

            sDir2.domain    = @(coor) isTop(coor);
            sDir2.direction = [2];
            sDir2.value     = perc*1;
            dir2 =  DirichletCondition(obj.mesh, sDir2);

            sDir3.domain    = @(coor) isTop(coor);
            sDir3.direction = [1];
            sDir3.value     = 0;
            dir3 =  DirichletCondition(obj.mesh, sDir3);

            s.dirichletFun = [dir, dir2, dir3];
%             s.dirichletFun = [dir, dir2];
            s.pointloadFun = [];%DistributedLoad(obj.mesh, sPL);

            Fext = zeros(obj.mesh.nnodes,2);

            obj.FextInitial = Fext; 
            
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

        function rshp = reshapeToVector(obj, A)
            rshp = reshape(A',[obj.uFun.nDofs,1]);
        end

        function rshp = reshapeToMatrix(obj, A)
            rshp = reshape(A,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

    end

end
