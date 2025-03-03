classdef levelSetTesting2 < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        nSubdomains
        tolNode
        Kel
        U
        
        mesh
        discMesh
        fineDiscMesh

        meshDomain
        boundaryConditions
        bcApplier

    end

    methods
        function obj = levelSetTesting2()
            clc; close all;

            obj.init();
            
            obj.createMesh();

            obj.coarseElasticProblem();

        end


        function init(obj)
            obj.nSubdomains = [3,1];
            obj.tolNode     = 1e-14;

            fileName = "UL_r0_1-P1.mat";
            data     = load(fileName);
            obj.Kel  = data.L;
            obj.U  = data.U;
        end

    end

    methods (Access = private)

        function createMesh(obj)
            mR  = obj.createStructuredMesh();        %creates reference fine mesh
            mRC = obj.createReferenceCoarseMesh(mR); %creates reference coarse mesh
            
            [mD,mSb,iC,lG,iCR, disMesh] = obj.createMeshDomain(mRC);
            [~,~,~,~,~,obj.fineDiscMesh] = obj.createMeshDomain(mR);
            
            % mD       -> simplified joint mesh
            % mSb      -> each element of the mesh domain in a cell
            % discMesh -> each domain of the mesh together but not joint

            obj.mesh     = mD;
            obj.discMesh = disMesh;
        end

        function mS = createStructuredMesh(obj)
             %UnitMesh better
            x1       = linspace(-1,1,50);
            x2       = linspace(-1,1,50);
            [xv,yv]  = meshgrid(x1,x2);
            [F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            bgMesh   = Mesh.create(s);

            lvSet    = obj.createLevelSetFunction(bgMesh);
            uMesh    = obj.computeUnfittedMesh(bgMesh,lvSet);
            mS       = uMesh.createInnerMesh();

        end
        
        function levelSet = createLevelSetFunction(~,bgMesh)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0;
            sLS.yCoorCenter = 0;
            sLS.radius      = 0.1;
            g               = GeometricalFunction(sLS);
            lsFun           = g.computeLevelSetFunction(bgMesh);
            levelSet        = lsFun.fValues;

        end

        function uMesh = computeUnfittedMesh(~,bgMesh,levelSet)
            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(levelSet);

        end

        function cMesh = createReferenceCoarseMesh(~,mR)
            xmax       = max(mR.coord(:,1));
            xmin       = min(mR.coord(:,1));
            ymax       = max(mR.coord(:,2));
            ymin       = min(mR.coord(:,2));
            coord(1,1) = xmin;
            coord(1,2) = ymin;
            coord(2,1) = xmax;
            coord(2,2) = ymin;
            coord(3,1) = xmax;
            coord(3,2) = ymax;
            coord(4,1) = xmin;
            coord(4,2) = ymax;
            
            connec   = [1 2 3 4];
            s.coord  = coord;
            s.connec = connec;
            cMesh    = Mesh.create(s);

        end

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode   = obj.tolNode;
            m               = MeshCreatorFromRVE(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end

        % --------------------

        function coarseElasticProblem(obj)
            dispFun                = LagrangianFunction.create(obj.mesh, obj.mesh.ndim,'P1');
            obj.boundaryConditions = obj.createCoarseBoundaryConditions();
            obj.createBCapplier();

            [lhs,LHSr]  = obj.computeStiffnessMatrix(dispFun);
            RHS         = obj.computeForces(lhs, dispFun);
            uRed        = LHSr\RHS;
            uCoarse     = obj.bcApplier.reducedToFullVectorDirichlet(uRed);

            u = obj.reconstructSolution(uCoarse, dispFun);
            
            EIFEMtesting.plotSolution(uCoarse, obj.mesh, 2, 3, 0, [], 0);     %plot solution with coarse mesh
            EIFEMtesting.plotSolution(u, obj.fineDiscMesh, 2, 4, 0, [], 0);     %plot solution with fine mesh

        end

        function bc = createCoarseBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            tol     = 1e-10;

            left    = @(coor) abs(coor(:,1)-xMin) <= tol;
            right   = @(coor) abs(coor(:,1)-xMax) <= tol;

            sDir{1}.domain    = @(coor) left(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;
           
            sPL{1}.domain     = @(coor) right(coor);
            sPL{1}.direction  = 2;
            sPL{1}.value      = -1;

            dirichletFun = [];

            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = cat(2, dirichletFun, dir);
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];

            for i = 1:numel(sPL)
                pl           = PointLoad(obj.mesh, sPL{i});
                pointloadFun = cat(2, pointloadFun, pl);

            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc             = BoundaryConditions(s);

        end

        function createBCapplier(obj)
            s.mesh                  = obj.meshDomain;
            s.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(s);

        end

        function [LHS,LHSr] = computeStiffnessMatrix(obj, dispFun)
            K = repmat(obj.Kel,[1,1, obj.mesh.nelem]);

            trial     = dispFun;
            test      = trial;
            s.fun     = dispFun;
            assembler = AssemblerFun(s);
            LHS       = assembler.assemble(K, test, trial);
            LHSr      = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);

        end

        function RHS = computeForces(obj,stiffness,u)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
            s.dim.ndofs = u.nDofs;
            s.BC        = obj.boundaryConditions;
            s.mesh      = obj.mesh;
            RHSint      = RHSintegrator.create(s);
            rhs         = RHSint.compute();
            % Perhaps move it inside RHSint?
            R           = RHSint.computeReactions(stiffness);
            RHS         = rhs+R;
            RHS         = obj.bcApplier.fullToReducedVectorDirichlet(RHS);

        end

        function u = reconstructSolution(obj, uCoarse, dispFun)
            nElem    = obj.mesh.nelem;
            dofConec = dispFun.getDofConnec();

            for ielem = 1:nElem
                uCelem     = uCoarse(dofConec(ielem,:));
                u(:,ielem) = obj.U*uCelem;
            end

            u = [u(:,1); u(:,2)];

        end


    end



end