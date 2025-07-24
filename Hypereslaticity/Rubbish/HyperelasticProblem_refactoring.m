classdef HyperelasticProblem_refactoring < handle

    properties (Access = public)
        uFun
        rFun
    end

    properties (Access = private)
        mesh
        neohookeanFun, linearElasticityFun
        material, materialElastic
        bc_case
        fileName
        meshGen
        printing
        bcApplier
        freeDofs, constrainedDofs
        refMesh
        nSubdomains
        tolSameNode
        boundaryConditions
        EIFEMp
%         Meifem
    end

    methods (Access = public)

        function obj = HyperelasticProblem_refactoring(cParams)
            close all;
            obj.init(cParams);
            tol = 1e-8;

            u  = obj.uFun;
            [hess, ~] = obj.computeHessian(u);
            hessf = @(x) hess*x;
            Meifem = @(r) obj.EIFEMp.apply(r);
            Milu      = obj.createILUpreconditioner(hess);
            Mmult     = @(r) Preconditioner.multiplePrec(r,Milu,Meifem,Milu,hessf);
%             coarseBasis = mat2cell(obj.EIFEMp.getProjectionMatrix(),[]);
            coarseBasis = num2cell(obj.EIFEMp.getProjectionMatrix(), 1);
            functionType = {'P1','P1','P1','P1','P1','P1','P1','P1'};
            Tfun =  ModalFunction.create(obj.refMesh,coarseBasis,functionType);
            f = animatedline;

            nsteps = cParams.nsteps;
            iter = 1;
            nIterPerStep = [];
            Usol = zeros(size(hess,1),1);
            reactions = zeros(u.nDofs,1);
            for iStep = 1:nsteps
                loadPercent = iStep/nsteps;
                [bc,~] = obj.createBoundaryConditions(loadPercent);
                % u  = obj.computeInitialElasticGuess(bc);
                u  = obj.applyDirichletToUFun(u,bc);

                intEnergy(iter) = obj.neohookeanFun.compute(u);
                %intEnergy(iter) = obj.linearElasticityFun.compute(u);
                Res = obj.computeResidual(u,loadPercent);
                resi(iter) = norm(Res);
                cgiters(iStep)=0;
                hasNotConverged = 1;
                while hasNotConverged
                    uVal = obj.reshapeToVector(u.fValues);
                    [hess, KR] = obj.computeHessian(u);
                    hessf = @(x) hess*x;
                    Res  = obj.computeResidual(u,loadPercent);
%                     incU = hess\(-Res);
                    [incU,residualPCG,errPCG,errAnormPCG] = PCG.solve(hessf,-Res,uVal(obj.freeDofs),Mmult,tol,Usol,obj.mesh,obj.bcApplier);
                    cgiters(iStep) =  cgiters(iStep)+length(residualPCG);
                    uVal(obj.freeDofs) = uVal(obj.freeDofs) + incU;
                    u.setFValues(obj.reshapeToMatrix(uVal));
                    reacFun = obj.computeReactions(KR,u);

                    iter = iter+1;
                    intEnergy(iter) = obj.neohookeanFun.compute(u);
                    Res        = obj.computeResidual(u,loadPercent);
                    resi(iter) = norm(Res);
                    hasNotConverged = resi(iter) > tol;
%                     hasNotConverged = obj.hasNotConverged(intEnergy);
                end
                normDisp(iStep) = Norm(u,'L2');
                normReac(iStep) = Norm(reacFun,'L2');
                obj.printFile(iStep, u, reacFun);              
                nIterPerStep(iStep) = obj.computeNumberOfIterations(iter,nIterPerStep,iStep);
                avgCGiters = cgiters./nIterPerStep;
                obj.plotStep(intEnergy,nIterPerStep,iStep, normDisp,normReac,avgCGiters);
                obj.rFun = reacFun;
            end

        end

    end

    methods (Access = private)

        function printFile(obj,iStep, u, reac)
            if obj.printing
                fun = {u, reac};
                funNames = {'Displacement', 'Reactions'};
                a.mesh     = obj.mesh;
                a.filename = ['SIM_',obj.fileName,'_',int2str(iStep)];
                a.fun      = fun;
                a.funNames = funNames;
                a.type     = 'Paraview';
                pst = FunctionPrinter.create(a);
                pst.print();
            end
        end
        
        function reacFun = computeReactions(obj,KR,u)
            uVal = obj.reshapeToVector(u.fValues);
            reactions = zeros(size(uVal));
            reactions(obj.constrainedDofs,1) = -KR*uVal;
            s.fValues = obj.reshapeToMatrix(reactions);
            s.order = 'P1';
            s.mesh  = obj.mesh;
            reacFun = LagrangianFunction(s);
        end

        function plotStep(obj,intEnergy,nIterPerStep,iStep, normDisp, normReac,cgiters)

            f = figure(1);
            clf(f)
            subplot(1,4,1)
            plot(intEnergy)
            title('Energy')
            xlabel('iteration')
            ylabel('Energy')
            subplot(1,4,2)
            plot(normDisp, normReac)
            title('Displacement-reaction')
            xlabel('displacement')
            ylabel('reactions')
            hold on
            subplot(1,4,3)
            bar(1:iStep, nIterPerStep)
%             title('Number of iterations to converge $\Delta$ F')
             title('Newton iterations')
            xlabel('step')
            ylabel('num. iterations per step')
            subplot(1,4,4)
            bar(1:iStep, cgiters)
            title('CG iterations')
            xlabel('step')
            ylabel('num. iterations per step')
            hold on
            drawnow
        end

        function hasNot = hasNotConverged(obj,energy)
            TOL = 10e-12;
            energyOld = energy(end-1);
            energy    = energy(end);
            hasNot = abs(energy-energyOld)/energy > TOL;
        end

        function init(obj,cParams)
            obj.printing    = cParams.printing;
            obj.bc_case     = cParams.bcCase;
            obj.fileName    = cParams.fileName;
            obj.meshGen     = cParams.meshGen;
            obj.nSubdomains = cParams.nSubdomains;
            obj.tolSameNode = 1e-10;
            obj.createReferenceMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMesh(obj.refMesh);
            obj.mesh = mD;            
            obj.createMaterial();
            obj.createDisplacementFun();
            obj.createFunctionals();
            obj.computeFreeDofs();
            [obj.boundaryConditions,dir] = obj.createBoundaryConditions(1);
            obj.createBCapplier();
            bS  = obj.refMesh.createBoundaryMesh();
            obj.createEIFEMPreconditioner(obj.refMesh,dir,iC,lG,bS,iCR,discMesh);
        end

        function u = computeInitialElasticGuess(obj, bc)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.materialElastic;
            s.dim = '2D';
            s.boundaryConditions = bc;
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            fem.solve();
            %             u = obj.reshapeToVector(fem.uFun.fValues);
            u = fem.uFun;
        end

        function createFunctionals(obj)
            obj.createLinearElasticityFunctional();
            obj.createNeohookeanFunctional();
        end

        function createLinearElasticityFunctional(obj)
            s.material = obj.materialElastic;
            s.mesh     = obj.mesh;
            elas = LinearElasticityFunctional(s);
            obj.linearElasticityFun = elas;
        end

        function createNeohookeanFunctional(obj)
            s.material = obj.material;
            s.mesh     = obj.mesh;
            neo = NeohookeanFunctional(s);
            obj.neohookeanFun = neo;
        end

        function [Hff,Hr] = computeHessian(obj,u)
            H = obj.neohookeanFun.computeHessian(u);
            %H = obj.linearElasticityFun.computeHessian(u);
            f = obj.freeDofs;
            r = obj.constrainedDofs;
            Hff = H(f,f);
            Hr  = H(r,:);
        end

        function Rf = computeResidual(obj,u,perc)
            Fext          = obj.computeExternalForces(perc);
            [Fint,FintEl] = obj.computeInternalForces(u);
            R = Fint - Fext;% - react;
            f = obj.freeDofs;
            Rf = R(f);
        end

        function nIter = computeNumberOfIterations(obj,iter,iterOld,iStep)
            if iStep == 1
                nIter = iter-1;
            else
                nIter = iter-1-sum(iterOld(1:(iStep-1)));
            end
        end

        function computeFreeDofs(obj)
            bc = obj.createBoundaryConditions(1);
            nDofs = obj.uFun.nDofs;
            obj.freeDofs = setdiff(1:nDofs,bc.dirichlet_dofs);
            obj.constrainedDofs = bc.dirichlet_dofs;
        end

        function createReferenceMesh(obj)
            switch obj.meshGen
                case {'Hole', 'HoleDirich'}
                    IM = Mesh.createFromGiD('hole_mesh_quad.m');
                    obj.refMesh = IM;
                    s.coord = obj.refMesh.coord;
                    s.connec = obj.refMesh.connec;
                     maxC = max(s.coord);
                    minC = min(s.coord);
                    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:) =...
                    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:)-1e-8;

                    s.coord(abs(s.coord(:,1)- minC(1))<1e-5 & abs(s.coord(:,2)- maxC(2))<1e-5,:) =...
                    s.coord(abs(s.coord(:,1)- minC(1))<1e-5 & abs(s.coord(:,2)- maxC(2))<1e-5,:)+[1e-8,-1e-8];

                    s.coord(abs(s.coord(:,1)- minC(1))<1e-5 & abs(s.coord(:,2)- minC(2))<1e-5,:) =...
                    s.coord(abs(s.coord(:,1)- minC(1))<1e-5 & abs(s.coord(:,2)- minC(2))<1e-5,:)+[1e-8,1e-8];
                    
                    s.coord(abs(s.coord(:,1)- maxC(1))<1e-5 & abs(s.coord(:,2)- minC(2))<1e-5,:) =...
                    s.coord(abs(s.coord(:,1)- maxC(1))<1e-5 & abs(s.coord(:,2)- minC(2))<1e-5,:)+[-1e-8,1e-8];
                    %                     s.coord(min(s.coord)) = s.coord(max(s.coord))-1-e9;
                    obj.refMesh = Mesh.create(s);
                case {'Bending', 'Traction'}
                    obj.refMesh = UnitQuadMesh(20,20);
                case {'Metamaterial'}
                    load('NegPoissMesh.mat')
                    s.coord = NegPoissMesh.coord;
                    s.connec = NegPoissMesh.connec;
                    maxC = max(s.coord);
                    minC = min(s.coord);
%                     s.coord(abs(s.coord(:,1)- maxC(1))<1e-5 & abs(s.coord(:,2)- maxC(2))<1e-5,:) =...
%                     s.coord(abs(s.coord(:,1)- maxC(1))<1e-5 & abs(s.coord(:,2)- maxC(2))<1e-5,:)-1e-3;

                    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:) =...
                    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:)-1e-2;

                    s.coord(abs(s.coord(:,1)- minC(1))<1e-5 & abs(s.coord(:,2)- maxC(2))<1e-5,:) =...
                    s.coord(abs(s.coord(:,1)- minC(1))<1e-5 & abs(s.coord(:,2)- maxC(2))<1e-5,:)+[1e-2,-1e-2];

                    s.coord(abs(s.coord(:,1)- minC(1))<1e-5 & abs(s.coord(:,2)- minC(2))<1e-5,:) =...
                    s.coord(abs(s.coord(:,1)- minC(1))<1e-5 & abs(s.coord(:,2)- minC(2))<1e-5,:)+[1e-2,1e-2];
                    
                    s.coord(abs(s.coord(:,1)- maxC(1))<1e-5 & abs(s.coord(:,2)- minC(2))<1e-5,:) =...
                    s.coord(abs(s.coord(:,1)- maxC(1))<1e-5 & abs(s.coord(:,2)- minC(2))<1e-5,:)+[-1e-2,1e-2];
                    %                     s.coord(min(s.coord)) = s.coord(max(s.coord))-1-e9;
                    obj.refMesh = Mesh.create(s);
                case 'EIFEMMesh'
                    load(obj.fileName)
                    s.coord  = EIFEoper.MESH.COOR;          
                    s.connec = EIFEoper.MESH.CN;
                    mS       = Mesh.create(s);
                    obj.refMesh = mS;
                otherwise
                    obj.refMesh = HexaMesh(2,1,1,20,5,5);
                    % obj.mesh = UnitHexaMesh(15,15,15);
            end
        end


        function [mD,mSb,iC,lG,iCR,discMesh] = createMesh(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end

        function createMaterial(obj)
            obj.material.mu = 1;
            obj.material.lambda = 1*1;
            obj.materialElastic = obj.createElasticMaterial();
        end

        function mat = createElasticMaterial(obj)
            G = obj.material.mu;
            L = obj.material.lambda;
            N = obj.mesh.ndim;
            K = 2/N*G + L;
            E1  = Isotropic2dElasticMaterial.computeYoungFromShearAndBulk(G,K,N)  ;
            nu1 = Isotropic2dElasticMaterial.computePoissonFromFromShearAndBulk(G,K,N);
            E2        = G*(3*L+2*G)/(L+G);
            nu2       = L / (2*(L+G));
            E         = ConstantFunction.create(E1,obj.mesh);
            nu        = ConstantFunction.create(nu1,obj.mesh);
            s.pdim    = obj.mesh.ndim;
            s.nelem   = obj.mesh.nelem;
            s.mesh    = obj.mesh;
            s.young   = E;
            s.poisson = nu;
            s.type = 'ISOTROPIC';
            s.ndim = obj.mesh.ndim;
            mat = Material.create(s);
        end

        function [bcs,dir] = createBoundaryConditions(obj,perc)
            bc = HyperelasticityTestBCs(obj.bc_case, obj.mesh, perc);
            bcs = bc.boundaryConditions;
            dir = bc.dir;
        end

        function createDisplacementFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            obj.uFun.setFValues(obj.uFun.fValues + 0);
        end

        function u = applyDirichletToUFun(obj,u,bc)
            u_k = reshape(u.fValues',[u.nDofs,1]);
            u_k(bc.dirichlet_dofs) = bc.dirichlet_vals;
            u.setFValues(reshape(u_k,[obj.mesh.ndim,obj.mesh.nnodes])');
        end

        function [Fext] = computeExternalForces(obj,perc)
                        Fext = perc*obj.boundaryConditions.pointloadFun.fValues;
                        Fext = obj.reshapeToVector(Fext);
%             Fext = zeros(obj.uFun.nDofs,1);
        end

        function [intFor,intForel] = computeInternalForces(obj,u)
            intFor = obj.neohookeanFun.computeGradient(u);
            %intForel = obj.linearElasticityFun.computeGradient(u);
            intForel = 0;
        end

        function rshp = reshapeToVector(obj, A)
            rshp = reshape(A',[obj.uFun.nDofs,1]);
        end

        function rshp = reshapeToMatrix(obj, A)
            rshp = reshape(A,[obj.mesh.ndim,obj.mesh.nnodes])';
        end

        function plotVector(obj, vect)
            s.fValues = obj.reshapeToMatrix(vect);
            s.mesh = obj.mesh;
            s.order = 'P1';
            a = LagrangianFunction(s);
            a.plot;
        end


        % Other useful stuff later on
        function lambdas = computeStretches(obj)
            I33 = Identity(obj.uFun);
            F = I33 + Grad(obj.uFun)';
            C = F'*F;
            lambdas = (Eigen(C).^0.5);
        end

        function sigma = computeCauchyStress(obj)
            mu = obj.material.mu;
            lambda = obj.material.lambda;
            I33 = Identity(obj.uFun);
            F = I33 + Grad(obj.uFun)';
            b = F*F';
            jac = Det(F);
            sigma = mu.*(b-Id)./jac + lambda.*(log(jac)).*Id./jac;
            sigma.ndimf = [2 2];
        end

        function  createEIFEMPreconditioner(obj,mR,dir,iC,lG,bS,iCR,dMesh)
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            EIFEMfilename = obj.fileName;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            filename        = EIFEMfilename;
            s.RVE           = TrainedRVE(filename);
            s.mesh          = obj.createCoarseMesh(mR);
%            s.mesh          = obj.loadCoarseMesh(mR);
            s.DirCond       = dir;
            s.nSubdomains = obj.nSubdomains;
            eifem           = EIFEM(s);


            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver = eifem;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'EIFEM';
            obj.EIFEMp = Preconditioner.create(ss);
%             Meifem = @(r) obj.EIFEMp.apply(r);
        end

        function mCoarse = createCoarseMesh(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh(mR);
%             s.meshReference = obj.loadReferenceCoarseMesh(mR);
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE(s);
            [mCoarse,~,~] = mRVECoarse.create();
        end

        function cMesh = createReferenceCoarseMesh(obj,mR)
            xmax = max(mR.coord(:,1));
            xmin = min(mR.coord(:,1));
            ymax = max(mR.coord(:,2));
            ymin = min(mR.coord(:,2));
            coord(1,1) = xmin;
            coord(1,2) = ymin;
            coord(2,1) = xmax;
            coord(2,2) = ymin;
            coord(3,1) = xmax;
            coord(3,2) = ymax;
            coord(4,1) = xmin;
            coord(4,2) = ymax;
            %             coord(1,1) = xmax;
            %             coord(1,2) = ymin;
            %             coord(2,1) = xmax;
            %             coord(2,2) = ymax;
            %             coord(3,1) = xmin;
            %             coord(3,2) = ymax;
            %             coord(4,1) = xmin;
            %             coord(4,2) = ymin;
%             connec = [1 2 3 4];
            connec = [2 3 4 1];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end

        function Milu = createILUpreconditioner(obj,LHS)
            s.LHS = LHS;
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);
        end

         function createBCapplier(obj)
            s.mesh                  = obj.mesh;
            s.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(s);
         end

          function d = createDomainDecompositionDofManager(obj,iC,lG,bS,mR,iCR)
            s.nSubdomains     = obj.nSubdomains;
            s.interfaceConnec = iC;
            s.interfaceConnecReshaped = iCR;
            s.locGlobConnec   = lG;
            s.nBoundaryNodes  = bS{1}.mesh.nnodes;
            s.nReferenceNodes = mR.nnodes;
            s.nNodes          = obj.mesh.nnodes;
            s.nDimf           = obj.mesh.ndim;
            d = DomainDecompositionDofManager(s);
        end

    end

end
