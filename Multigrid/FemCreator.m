classdef FemCreator < handle

    properties (Access = public)
        LHS
        RHS
        bc
        bcApplier
        solver
    end

    properties (Access = private)
        coarseMeshes
        nDimf
        nLevel
        coarseBc
        coarseDispFun
        coarseMaterial
        type
        scale
        pdim
        tol
        solverCase
    end

    methods (Access = public)
        function obj = FemCreator(cParams)
            obj.init(cParams);
            u             = obj.createLevelLagrangianFunction();
            bc            = obj.createLevelBoundaryConditions();
            obj.bcApplier = obj.createLevelBoundaryConditionsApplier(bc);
            mat           = obj.createLevelMaterial();

            LHSf              = obj.createLHS(u,bc,mat);
            RHSf              = obj.createRHS(u,bc,mat,LHSf);      
            [obj.LHS,obj.RHS] = obj.computeReducedSystem(LHSf,RHSf);
 
            obj.createSolver();
            obj.createSolverCoarse();
        end
    end

    methods (Access = private)
        function init (obj,cParams)
            obj.coarseMeshes = cParams.coarseMeshes;
            obj.nDimf        = cParams.nDimf;
            obj.nLevel       = cParams.nLevel;
            obj.type         = cParams.type;
            obj.scale        = cParams.scale;
            obj.pdim         = cParams.pdim;
            obj.tol          = cParams.tol;
            obj.solverCase   = cParams.solverCase;
        end

        function u = createLevelLagrangianFunction(obj)
            for i = 1:obj.nLevel
                u(i) = LagrangianFunction.create(obj.coarseMeshes{i}, obj.nDimf,'P1');
            end
        end

        function bcApplier = createLevelBoundaryConditionsApplier(obj,bc)
            for i = 1:obj.nLevel
                s.mesh                 = obj.coarseMeshes{i};
                s.boundaryConditions   = bc(i);
                bcApplier(i)           = BCApplier(s);
            end
        end

        function bc = createLevelBoundaryConditions(obj)
            [Dir,PL] = obj.createRawBoundaryConditions();
            for i = 1:obj.nLevel
                bc(i) = obj.createBoundaryConditions(obj.coarseMeshes{i},Dir,PL);
            end
        end

        function bc = createBoundaryConditions(obj,mesh,Dir,PL)
            dirichlet = DirichletCondition(mesh,Dir);
            pointload = PointLoad(mesh,PL);
             % need this because force applied in the face not in a point
            pointload.values=pointload.values/size(pointload.dofs,1);
            fvalues = zeros(mesh.nnodes*obj.nDimf,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues = reshape(fvalues,obj.nDimf,[])';
            pointload.fun.setFValues(fvalues);

            s.pointloadFun = pointload;
            s.dirichletFun = dirichlet;
            s.periodicFun =[];
            s.mesh = mesh;
            bc          = BoundaryConditions(s);
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir.domain    = @(coor) isLeft(coor);
            Dir.direction = [1,2,3];
            Dir.value     = 0;

            PL.domain    = @(coor) isRight(coor);
            PL.direction = 2;
            PL.value     = -1;
        end

        function mat = createLevelMaterial(obj)
            for i = 1:obj.nLevel
                mat(i) = obj.createMaterial(obj.coarseMeshes{i});
            end
        end

        function [young,poisson] = computeElasticProperties(obj,mesh)
            E  = 1;
            nu = 1/3;
            young   = ConstantFunction.create(E,mesh);
            poisson = ConstantFunction.create(nu,mesh);      
        end

        function material = createMaterial(obj,mesh)
            [young,poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            tensor    = Material.create(s);
            material  = tensor;
        end

        function LHS = createLHS(obj,dispFunLevels,bcLevels,matLevels)
            for i = 1:obj.nLevel
                disFun     = dispFunLevels(i);
                mat        = matLevels(i);
                m          = obj.coarseMeshes{i};
                LHS{i}     = obj.computeStiffnessMatrix(m,mat,disFun);
            end
        end

        function RHS = createRHS(obj,dispFunLevels,bcLevels,matLevels,LHSlevels)
            for i = 1:obj.nLevel
                dispFun    = dispFunLevels(i);
                bc         = bcLevels(i);
                mat        = matLevels(i);
                LHS        = LHSlevels{i};
                m          = obj.coarseMeshes{i};
                RHS{i}     = obj.computeForces(m,dispFun,bc,material,LHS,obj.solverCase);
            end
        end

        function Kred = computeKred(obj,m,mat,u,bc)
            K    = obj.computeStiffnessMatrix(m,mat,u);
            Kred = bc.fullToReducedMatrix(K);
        end
        
        function LHS = computeStiffnessMatrix(obj,mesh,material,displacementFun)
            ndimf      = displacementFun.ndimf;
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.test     = LagrangianFunction.create(mesh,ndimf, 'P1');
            s.trial    = displacementFun;
            s.material = material;
            s.quadratureOrder = 2;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function Fred = computeFred(obj,m,u,bc)
            b    = obj.computeRHS(m,u,bc);
            Fred = bc.fullToReducedVector(b);
        end

         function forces = computeForces(obj,mesh,dispFun,boundaryConditions,material,stiffness,solverCase)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
            s.dim.ndofs = dispFun.nDofs;
            s.BC       = boundaryConditions;
            s.mesh     = mesh;
            s.material = material;
%             s.globalConnec = mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            % Perhaps move it inside RHSint?
            if strcmp(solverCase,'REDUCED')
                R = RHSint.computeReactions(stiffness);
                forces = rhs+R;
            else
                forces = rhs;
            end
         end

         function [LHSred,RHSred] = computeReducedSystem(obj,LHSfull,RHSfull)
             for i = 1:obj.nLevel
                 LHSred{i} = obj.bcApplier(i).fullToReducedMatrixDirichlet(LHSfull{i});
                 RHSred{i} = obj.bcApplier(i).fullToReducedVectorDirichlet(RHSfull{i});
             end
         end

        function createSolver(obj)
            for i = 2:obj.nLevel+1
                s.maxIter             = 100;
                s.tol                 = obj.tol;
                s.solverType          = 'ITERATIVE';
                s.iterativeSolverType = 'CG';
                obj.solver{i}         = Solver.create(s);
            end
        end

        function createSolverCoarse(obj)
            s.maxIter             = 100000;
            s.tol                 = obj.tol;
            s.solverType          = 'ITERATIVE';
            s.iterativeSolverType = 'CG';

            obj.solver{1}         = Solver.create(s);
        end

    end
end

