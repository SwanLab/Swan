classdef levelSetTesting2 < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        nSubdomains
        fileName
        tol
        meshDomain
        bc
        bcApplier
        LHS
        mesh
        Kel
        U

        boundaryMesh
        boundaryMeshJoined
        localGlobalConnecBd
        dLambda
        displacementFun
    end

    methods
        function obj = levelSetTesting2()
            clc; close all;

            obj.init();
            
            mR  = obj.createStructuredMesh();
            mRc = obj.createReferenceCoarseMesh(mR);

            obj.coarseElasticProblemHere(mRc);

        end

        function init(obj)
            obj.nSubdomains = [2, 1];
            obj.fileName    = "UL_r0_3.mat";
            obj.tol         = 1e-14;
            data            = load(obj.fileName);
            obj.Kel         = data.L;
            obj.U           = data.U;

        end

        function mS = createStructuredMesh(obj)
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
            obj.mesh = mS;

            obj.boundaryMesh = obj.mesh.createBoundaryMesh();

            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.mesh.createSingleBoundaryMesh();

        end

        function levelSet = createLevelSetFunction(~,bgMesh)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0.5;
            sLS.yCoorCenter = 0.5;
            sLS.radius      = 0.3;
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

        function cMesh = createReferenceCoarseMesh(~, mR)
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
            
            connec   = [2 3 4 1];
            s.coord  = coord;
            s.connec = connec;
            cMesh    = Mesh.create(s);

        end

        function coarseElasticProblemHere(obj, mR)
            bS                     = mR.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR]     = obj.createMeshDomain(mR);
            obj.mesh               = mD;
            dispFun                = LagrangianFunction.create(obj.mesh, obj.mesh.ndim,'P1');
            obj.bc = createCoarseBoundaryConditions(obj);
            obj.createBCapplier();

            [LHS,LHSr] = obj.computeStiffnessMatrix(dispFun);
            RHS        = obj.computeForces(LHS, dispFun);
            c          = obj.computeCmat();
            [u, L]     = obj.computeDisplacementHere(c, LHS, RHS);


        end

        function [mD,mSb,iC,lG,iCR] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode   = obj.tol;
            m               = MeshCreatorFromRVE(s);

            [mD,mSb,iC,~,lG,iCR] = m.create();

        end

        function createBCapplier(obj)
            s.mesh                  = obj.meshDomain;
            s.boundaryConditions    = obj.bc;
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
            s.BC        = obj.bc;
            s.mesh      = obj.mesh;
            RHSint      = RHSintegrator.create(s);
            rhs         = RHSint.compute();
            R           = RHSint.computeReactions(stiffness);

            RHS = rhs+R;
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(RHS);
        end

        function Cg = computeCmat(obj)
            obj.createDisplacementFunHere();

            s.quadType = 2;
            s.mesh     = obj.boundaryMeshJoined;

            lhs    = LHSintegrator_ShapeFunction_fun(s);
            test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); % !!
            ndimf  = 2;
            Lx     = max(obj.mesh.coord(:,1)) - min(obj.mesh.coord(:,1));
            Ly     = max(obj.mesh.coord(:,2)) - min(obj.mesh.coord(:,2));

            f1 = @(x) [1/(4)*(1-x(1,:,:)).*(1-x(2,:,:));...
                    1/(4)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
            f2 = @(x) [1/(4)*(1+x(1,:,:)).*(1-x(2,:,:));...
                    1/(4)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
            f3 = @(x) [1/(4)*(1+x(1,:,:)).*(1+x(2,:,:));...
                    1/(4)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
            f4 = @(x) [1/(4)*(1-x(1,:,:)).*(1+x(2,:,:));...
                    1/(4)*(1-x(1,:,:)).*(1+x(2,:,:))  ];
            f  = {f1 f2 f3 f4}; %

            nfun = size(f,2);
            Cg   = [];

            for i = 1:nfun
                obj.dLambda{i}  = AnalyticalFunction.create(f{i},ndimf,obj.boundaryMeshJoined);
                Ce              = lhs.compute(obj.dLambda{i},test);
                [iLoc,jLoc,vals] = find(Ce);
    
                l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
                l2g_dof = l2g_dof(:);
                iGlob   = l2g_dof(iLoc);
                Cg      = [Cg sparse(iGlob,jLoc,vals, obj.displacementFun.nDofs, obj.dLambda{i}.ndimf)];
            end

        end

        function createDisplacementFunHere(obj)
            obj.displacementFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

        function [u, L] = computeDisplacementHere(obj,c,LHS,RHS)
            K   = LHS;
            nC  = size(c,2);
            Z   = zeros(nC);
            LHS = [K, c; c' Z];

            ud    = zeros(nC,1);
            ud(1) = 1;
            RHS   = [RHS; ud];

            sol = LHS\RHS;
            u   = sol(1:obj.displacementFun.nDofs);
            L   = sol(obj.displacementFun.nDofs+1:end); 

        end






    end



end