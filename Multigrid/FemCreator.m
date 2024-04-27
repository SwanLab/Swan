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
            pointload.values=pointload.values/size(pointload.dofs,1); % need this because force applied in the face not in a point
            s.pointloadFun = pointload;
            s.dirichletFun = dirichlet;
            s.periodicFun =[];
            s.mesh = mesh;
            bc          = BoundaryConditions(s);
            %             rawBc       = obj.createRawBoundaryConditions(mesh);
            %             dim         = getFunDims(obj,mesh,disp);
            %             rawBc.ndimf = dim.ndimf;
            %             rawBc.ndofs = dim.ndofs;
            %             s.mesh      = mesh;
            %             s.scale     = 'MACRO';
            %             s.bc        = {rawBc};
            %             s.ndofs     = dim.ndofs;
            %             bc          = BoundaryConditions(s);
            %             bc.compute();
        end

        function dim = getFunDims(obj,mesh,disp)
            d.ndimf     = disp.ndimf;
            d.nnodes    = size(disp.fValues, 1);
            d.ndofs     = d.nnodes*d.ndimf;
            d.nnodeElem = mesh.nnodeElem;
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim         = d;
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            %             dirichletNodes            = abs(mesh.coord(:,1)-0) < 1e-12;
            %             rightSide                 = max(mesh.coord(:,1));
            %             isInRight                 = abs(mesh.coord(:,1)-rightSide)< 1e-12;
            %             forceNodes                = isInRight;
            %             nodes                     = 1:mesh.nnodes;
            %             bcDir                     = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
            %             nodesdir                  = size(nodes(dirichletNodes),2);
            %             bcDir(1:nodesdir,end+1)   = 1;
            %             bcDir(nodesdir+1:end,end) = 2;
            %             bcDir(:,end+1)            = 0;
            %             bc.dirichlet              = bcDir;
            %             bc.pointload(:,1)         = nodes(forceNodes);
            %             bc.pointload(:,2)         = 2;
            %             bc.pointload(:,3)         = -1/length(forceNodes);

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
            E1  = 1;
            nu1 = 1/3;
            E   = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,mesh);
            nu  = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,mesh);
            young   = E;
            poisson = nu;
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

        %         function mat = createMaterial(obj,mesh)
        %             s.mesh  = mesh;
        %             s.type  = obj.type;
        %             s.scale = obj.scale;
        %             ngaus   = 1;
        %             Id      = ones(mesh.nelem,ngaus);
        %             s.ptype = 'ELASTIC';
        %             s.pdim  = obj.pdim;
        %             s.nelem = mesh.nelem;
        %             s.mesh  = mesh;
        %             s.kappa = .9107*Id;
        %             s.mu    = .3446*Id;
        %             mat     = Material.create(s);
        %             mat.compute(s);
        %         end

        function LHS = createLHS(obj,dispFunLevels,bcLevels,matLevels)
            for i = 1:obj.nLevel
                disFun     = dispFunLevels(i);
%                 bc         = bcLevels(i);
                mat        = matLevels(i);
                m          = obj.coarseMeshes{i};

                LHS{i}        = obj.computeStiffnessMatrix(m,mat,disFun);
                %                 obj.LHS    = obj.bcApplier(i).fullToReducedMatrixDirichlet(LHS);
                %                 obj.LHS{i} = obj.computeKred(m,mat,disFun,bc);
            end
        end

        function RHS = createRHS(obj,dispFunLevels,bcLevels,matLevels,LHSlevels)
            for i = 1:obj.nLevel
                dispFun     = dispFunLevels(i);
                bc         = bcLevels(i);
                mat        = matLevels(i);
                LHS        = LHSlevels{i};
                m          = obj.coarseMeshes{i};

                RHS{i}     = obj.computeForces(m,dispFun,bc,material,LHS,obj.solverCase);
%                 RHS{i} = obj.computeFred(m,u,bc,LHS);
            end
        end

        function Kred = computeKred(obj,m,mat,u,bc)
            K    = obj.computeStiffnessMatrix(m,mat,u);
            Kred = bc.fullToReducedMatrix(K);
        end
        %
        %         function k = computeStiffnessMatrix(obj,mesh,material,displacementFun)
        %             s.type     = 'ElasticStiffnessMatrix';
        %             s.mesh     = mesh;
        %             s.fun      = displacementFun;
        %             s.material = material;
        %             lhs        = LHSintegrator.create(s);
        %             k          = lhs.compute();
        %         end

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

%         function RHS = computeRHS(obj,mesh,dispFun,boundaryConditions)
%             dim.ndimf     = dispFun.ndimf;
%             dim.nnodes    = size(dispFun.fValues, 1);
%             dim.ndofs     = dim.nnodes*dim.ndimf;
%             dim.nnodeElem = mesh.nnodeElem;
%             dim.ndofsElem = dim.nnodeElem*dim.ndimf;
%             c.dim         = dim;
%             c.mesh        = mesh;
%             c.BC          = boundaryConditions;
%             RHS           = RHSintegrator_ElasticMacro(c);
%             RHS           = RHS.compute();
%         end


         function forces = computeForces(obj,mesh,dispFun,boundaryConditions,material,stiffness,solverCase)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
%             s.dim       = obj.getFunDims();
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
                s.maxIter             = 5;
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

