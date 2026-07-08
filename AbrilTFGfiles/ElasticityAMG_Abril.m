classdef ElasticityAMG_Abril < handle

% WARNING!
% - Install a valid python version (3.8 for example)
% - Install pyAMG python library
% - Install further python dependencies (numpy, scipy, ...)

    properties(Access=public)
        uSol
        residual
        error
        errAnorm
    end

    properties (Access = private)
        mesh
        young
        poisson
        material
        stateProblem
        nSubdomains
        boundaryConditions
        bcApplier
        fileName
    end

    methods (Access = public)

        function obj = ElasticityAMG_Abril(cParams)
            obj.init(cParams);
            obj.createMesh();
            obj.computeElasticProperties();
            obj.createMaterial();
            obj.solveElasticProblem();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nSubdomains = cParams.nSubdomains;
            obj.fileName    = cParams.fileNameEIFEM;

        end

        function createMesh(obj)
            load(obj.fileName,'mesh');
            s.coord    = mesh.coord;
            s.connec   = mesh.connec;

            xmin = min(s.coord(:,1));            
            xmax = max(s.coord(:,1));
            ymin = min(s.coord(:,2));
            ymax = max(s.coord(:,2));
            delta = 1e-8;
            s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymax,:) =...
                s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymax,:)+[-delta,-delta];
            s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymin,:) =...
                s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymin,:)+[-delta,delta];
            s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymax,:) =...
                s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymax,:)+[delta,-delta];
            s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymin,:) =...
                s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymin,:)+[delta,delta];
            mS = Mesh.create(s);

            s.nsubdomains   = obj.nSubdomains; %  num along x and y
            s.meshReference = mS;
            s.tolSameNode = 1e-9;
            m = MeshCreatorFromRVE.create(s);
            [mD,~,~,~,~,~,~] = m.create();
            obj.mesh      = mD;
        end

        function computeElasticProperties(obj)
            E  = 1;
            nu = 1/3;
            obj.young   = ConstantFunction.create(E,obj.mesh);
            obj.poisson = ConstantFunction.create(nu,obj.mesh);
        end

        function createMaterial(obj)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = obj.young;
            s.poisson = obj.poisson;
            tensor    = Material.create(s);
            obj.material = tensor;
        end

        function solveElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.material;
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            obj.boundaryConditions = s.boundaryConditions;
            obj.createBCApplier();
            
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = obj.createSolver(s);

            [LHS,RHS] = obj.createElasticProblem();
            Usol   = LHS\RHS;
            % Ufull  = obj.bcApplier.reducedToFullVectorDirichlet(Usol); 
            tic
            [obj.uSol,obj.residual,obj.error,obj.errAnorm] = s.solverCase.solve(LHS,RHS,Usol);
            % xFull = obj.bcApplier.reducedToFullVectorDirichlet(obj.uSol);
            t_AMG=toc

            % fem = ElasticProblem(s);
            % fem.solve();
            % obj.stateProblem = fem;
        end


        function createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bc = BCApplier(s);
            obj.bcApplier = bc;
        end

        function solver = createSolver(obj,s)
            BCAp = BCApplier(s);
            Rfull  = obj.computeRigidBodyModes([0.5,0.5]);
            for i = 1:size(Rfull,2)
                R(:,i) = BCAp.fullToReducedVectorDirichlet(Rfull(:,i));
            end
            s.type = 'ELASTIC';
            s.nullSpace = R;
            s.nLevels = 5;
            s.tol = 1e-8;
            s.maxIter = 1;
            p     = pyAMG.create(s);

            sS.preconditioner = p;
            sS.tol = 1e-8;
            solver = PCG_AMG(sS);
        end

        function R = computeRigidBodyModes(obj,refPoint)
            rigModes = RigidBodyFunction.create(obj.mesh,refPoint);
            RFun = rigModes.projectBasisFunctions('P1');
            for i = 1:length(RFun)
                R(:,i) = reshape(RFun{i}.fValues',[],1);
            end
        end

        function [LHSr,RHS,lhs]= createElasticProblem(obj)
            uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            [lhs,LHSr]= obj.computeStiffnessMatrix(uFun);
            RHS= obj.computeForces(lhs,uFun);
        end

        function [LHS,LHSr]= computeStiffnessMatrix(obj,uFun)
            C     = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            LHS = IntegrateLHS(f,uFun,uFun,obj.mesh,'Domain',2);
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);

        end

        function RHS =  computeForces(obj,stiffness,u)
            bc  = obj.boundaryConditions;
            t   = bc.tractionFun;
            rhs = zeros(u.nDofs,1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhsi = t(i).computeRHS(u);
                    rhs  = rhs + rhsi;
                end
            end
             bc      = obj.boundaryConditions;
             dirich  = bc.dirichlet_dofs;
             dirichV = bc.dirichlet_vals;
             if ~isempty(dirich)
                 R = -stiffness(:,dirich)*dirichV;
             else
                 R = zeros(sum(u.nDofs(:)),1);
             end
             rhs = rhs+R;
             RHS = obj.bcApplier.fullToReducedVectorDirichlet(rhs);
        end

        function bc = createBoundaryConditions(obj)
            minx = min(obj.mesh.coord(:,1));
            maxx = max(obj.mesh.coord(:,1));
            miny = min(obj.mesh.coord(:,2));
            maxy = max(obj.mesh.coord(:,2));

            tolBound = 1e-9;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);


            sDir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL.domain    = @(coor) isRight(coor);
            sPL.direction = [1];
            sPL.value     = [1];

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            % pointloadFun = [];
            % for i = 1:numel(sPL)
            %     pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
            %     pointloadFun = [pointloadFun, pl];
            % end
            pointload = TractionLoad(obj.mesh,sPL,'DIRAC');
            s.pointloadFun = pointload;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end

    end

end