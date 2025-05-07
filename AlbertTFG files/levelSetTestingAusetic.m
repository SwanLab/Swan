classdef levelSetTestingAusetic < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        nSubdomains
        fileNameEIFEM
        tolSameNode
        meshDomain
        boundaryConditions
        bcApplier
        LHS
        mesh
        Kel
        U
        RmeshFilename
    end

    methods
        function obj = levelSetTestingAusetic()
            clc; close all;
            obj.init();

            mR       = obj.createReferenceMesh();
            [mD,mSb,iC,lG,iCR, discMesh] = obj.createMeshDomain(mR);
            mRcoarse = obj.createReferenceCoarseMesh(mR);
            u = obj.createCoarseElasticProblem(mRcoarse);
            u = u(:);
            
            %CoarsePlotSolution(u, discMesh, [], "Ausetic Fine test")
            %EIFEMtesting.plotSolution(u, discMesh, 20, 1, 0, [], 0)
            

        end

        function init(obj)
            obj.nSubdomains  = [4 1]; %nx ny
            obj.fileNameEIFEM = "UL-Ausetic";
            obj.RmeshFilename = "DEF_Q4auxL_1.mat";
            obj.tolSameNode = 1e-14;
            data = load(obj.fileNameEIFEM);
            obj.Kel = data.L;
            obj.U = data.U;
        end
        
        function mS = createReferenceMesh(obj)
            mS = obj.loadAuseticMesh();
        end

        function mesh = loadAuseticMesh(obj)
            Data = load(obj.RmeshFilename);
            s.coord  = Data.EIFEoper.MESH.COOR;
            s.connec = Data.EIFEoper.MESH.CN;

            mesh = Mesh.create(s);

        end

        function mS = createStructuredMesh(obj)
             %UnitMesh better
            x1       = linspace(-1,1,20);
            x2       = linspace(-1,1,20);
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
        
        function [mD,mSb,iC,lG,iCR, discMesh] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE(s);
            [mD,mSb,iC,~,lG,iCR, discMesh] = m.create();
        end
        
        function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichletFun = [];
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.meshDomain, Dir{i});
                dirichletFun = cat(2,dirichletFun, dir);
            end

            pointload = PointLoad(mesh,PL);
            % need this because force applied in the face not in a point
            pointload.values        = pointload.values/size(pointload.dofs,1);
            fvalues                 = zeros(mesh.nnodes*mesh.ndim,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues                 = reshape(fvalues,mesh.ndim,[])';
            pointload.fun.fValues   = fvalues;

            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
        end
    
        function [Dir,PL] = createRawBoundaryConditions(obj)
            minx = min(obj.meshDomain.coord(:,1));
            maxx = max(obj.meshDomain.coord(:,1));
            %miny = min(obj.meshDomain.coord(:,2));
            maxy = max(obj.meshDomain.coord(:,2));
            tolBound = 1e-8;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            %isBottom = @(coor) (abs(coor(:,2) - miny)   < tolBound);
            isTop    = @(coor) (abs(coor(:,2) - maxy)   < tolBound);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir{1}.domain    = @(coor) isLeft(coor)| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

            %             Dir{2}.domain    = @(coor) isRight(coor) ;
            %             Dir{2}.direction = [2];
            %             Dir{2}.value     = 0;

            PL.domain    = @(coor) isTop(coor);
            PL.direction = 2;
            PL.value     = -0.1;
            %             PL.domain    = @(coor) isRight(coor);
            %             PL.direction = [1];
            %             PL.value     = [0.1];
        end
        
        function createBCapplier(obj)
            s.mesh                  = obj.meshDomain;
            s.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(s);
        end

        function [LHSr,RHSr,lhs] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix();
            RHSr       = obj.computeForces(lhs,u);
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

        function [young,poisson] = computeElasticProperties(~,mesh)
            E  = 1;
            nu = 1/3;
            young   = ConstantFunction.create(E,mesh);
            poisson = ConstantFunction.create(nu,mesh);
        end

        function [LHS,LHSr] = computeStiffnessMatrix(obj, dispFun)
            
            K = repmat(obj.Kel,[1,1, obj.mesh.nelem]);
            
            trial = dispFun;
            test  = trial;
            s.fun = dispFun;
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(K, test, trial);
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);
        end

        function RHS = computeForces(obj,stiffness,u)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
            s.dim.ndofs = u.nDofs;
            s.BC        = obj.boundaryConditions;
            s.mesh      = obj.mesh;
            RHSint      = RHSIntegrator.create(s);
            rhs         = RHSint.compute();
            % Perhaps move it inside RHSint?
            R           = RHSint.computeReactions(stiffness);
            RHS = rhs+R;
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(RHS);
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

            % xmax       = max(mR.coord(:,1));
            % xmin       = min(mR.coord(:,1));
            % ymax       = max(mR.coord(:,2));
            % ymin       = min(mR.coord(:,2));
            % coord(1,1) = xmin;
            % coord(1,2) = (ymin+ymax)/2;
            % coord(2,1) = (xmax+xmin)/2;
            % coord(2,2) = ymin;
            % coord(3,1) = xmax;
            % coord(3,2) = (ymin+ymax)/2;
            % coord(4,1) = (xmax+xmin)/2;
            % coord(4,2) = ymax;
            % 
            connec   = [1 2 3 4];
            s.coord  = coord;
            s.connec = connec;
            cMesh    = Mesh.create(s);
        end

        function [u] = createCoarseElasticProblem(obj, mRcoarse)
            bS                 = mRcoarse.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR, discMesh] = obj.createMeshDomain(mRcoarse);
            obj.mesh           = mD;
            

            dispFun                = LagrangianFunction.create(obj.mesh, obj.mesh.ndim,'P1');
            obj.boundaryConditions = createCoarseBoundaryConditions(obj);
            obj.createBCapplier();
            [LHS,LHSr]             = obj.computeStiffnessMatrix(dispFun);
            
            

            RHS = obj.computeForces(LHS, dispFun);
            
            uRed = LHSr\RHS;
            uCoarse = obj.bcApplier.reducedToFullVectorDirichlet(uRed);
            %CoarsePlotSolution(uCoarse, mD, [], "Ausetic Coarse test")
            %EIFEMtesting.plotSolution(uCoarse, mD, 1, 2, 0, [], 0)


            u = obj.reconstructSolution(uCoarse, dispFun);
            
            % for i = 1:obj.nSubdomains(1)
            %     EIFEMtesting.plotSolution(u(:,i), mSb{i}, 1, i, 0, [], 0)
            % 
            % 
            % end


            
        end

        function [lhs,RHSr,LHSr] = coarseElasticProblem(obj)
            
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial(obj.meshDomain);
            [lhs, LHSr] = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
            RHSr       = obj.computeForces(lhs,u);

        end

        function u = reconstructSolution(obj, uCoarse, dispFun)
            nElem = obj.mesh.nelem;
            dofConec = dispFun.getDofConnec();
            for ielem = 1:nElem
                uCelem = uCoarse(dofConec(ielem,:));
                u(:,ielem) =  obj.U*uCelem;
            end

        end
        
        function bc = createCoarseBoundaryConditions(obj)

            % v                    = zeros(8,1);
            % v(obj.nodeDirection) = 1;
            % nRes                 = [1 1 2 2 3 3 4 4]*v;
            % assignMatrix         = [2 1 0 0 0 0 0 0
            %                         0 0 2 1 0 0 0 0
            %                         0 0 0 0 2 1 0 0
            %                         0 0 0 0 0 0 2 1
            %                         1 2 1 2 1 2 1 2];
            % 
            % vSimp       = assignMatrix*v;
            % dirs        = cell(5,1);
            % [dirs{:,1}] = deal([1, 2]);
            % [dirs{nRes}]  = deal(vSimp(nRes));
            % [dirs{end}]   = deal(vSimp(end));

            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            tol     = 1e-10;

            left    = @(coor) abs(coor(:,1)-xMin) <= tol;
            right   = @(coor) abs(coor(:,1)-xMax) <= tol;
            % corner3 = @(coor) abs(coor(:,1)-xMax) <= tol & abs(coor(:,2)-yMax)<= tol;
            % corner4 = @(coor) abs(coor(:,1)-xMin) <= tol & abs(coor(:,2)-yMax)<= tol;

            sDir{1}.domain    = @(coor) left(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;
            % sDir{1}.ndim      = 3;
            
            % sDir{2}.domain    = @(coor) corner2(coor);
            % sDir{2}.direction = [1,2];
            % sDir{2}.value     = 0;
            % 
            % sDir{3}.domain    = @(coor) corner3(coor);
            % sDir{3}.direction = [1,2];
            % sDir{3}.value     = 0;
            % 
            % sDir{4}.domain    = @(coor) corner4(coor);
            % sDir{4}.direction = [1,2];
            % sDir{4}.value     = 0;
            
            % sDir{5}.domain    = @(coor) corner1(coor);
            % sDir{5}.direction = 1;
            % sDir{5}.value     = 1;
            

            % isForce = @(coor) (abs(coor(:,1) - xMin)   < 1e-10);


            sPL{1}.domain    = @(coor) right(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = 1;
            %sPL{1}.ndim      = 3;
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = cat(2, dirichletFun, dir);
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = cat(2, pointloadFun, pl);
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end

    end
end
