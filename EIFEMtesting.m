classdef EIFEMtesting < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        EIFEMfilename
        meshDomain
        boundaryConditions
        bcApplier
        solverCase
        material
        meshReference
        interfaceMeshReference
        meshSubDomain
        nSubdomains
        ninterfaces
        interfaceMeshSubDomain
        globalMeshConnecSubDomain
        interfaceMeshConnecSubDomain
        subDomainContact
        cornerNodes
        quad
        interfaceConnec
        locGlobConnec
        localGlobalDofConnec
        interfaceDof
        interfaceDom
        weight
        Lchol
        L
        U
        D

        eigenModes
        Kmodal
        MmodalPrecond
        displacementFun
        LHS
        RHS
        scale
        ndimf    
        functionType
        EIFEMsolver
        refLHS
        KeifemContinuous
        EIFEMprojection
    end

    methods (Access = public)

        function createFineMesh(obj)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.meshReference;         
            m = MeshCreatorFromRVE(s);
            [obj.meshDomain,obj.meshSubDomain,obj.interfaceConnec,~,obj.locGlobConnec] = m.create();
        end


        function obj = EIFEMtesting()
            close all
            obj.init();
            obj.createReferenceMesh();
            obj.createFineMesh();

            obj.displacementFun      = LagrangianFunction.create(obj.meshDomain, obj.ndimf,obj.functionType);
            [obj.boundaryConditions,Dir,PL] = obj.createBoundaryConditions(obj.meshDomain);
            ss.mesh                  = obj.meshDomain;
            ss.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier            = BCApplier(ss);
            obj.localGlobalDofConnec = obj.createlocalGlobalDofConnec();
            [obj.interfaceDof,obj.interfaceDom] = obj.computeLocalInterfaceDof();

%             obj.quad               = Quadrature.set(obj.meshDomain.type);
%             obj.quad.computeQuadrature('QUADRATIC');
%             obj.createDomainMaterial();
%               obj.computeForces();
    
            obj.material = obj.createMaterial(obj.meshDomain);
            obj.LHS      = obj.computeStiffnessMatrix(obj.meshDomain,obj.displacementFun,obj.material);
            obj.RHS      = obj.computeForces(obj.LHS);
           

            meshDomainCoarse = obj.createCoarseMesh();

            refDisp     = LagrangianFunction.create(obj.meshReference, obj.ndimf,obj.functionType);
            refMat      = obj.createMaterial(obj.meshReference);
            obj.refLHS      = obj.computeStiffnessMatrix(obj.meshReference,refDisp,refMat);

            obj.EIFEMsolver = obj.createEIFEM(meshDomainCoarse,Dir,obj.refLHS);
%             obj.computeKEIFEMglobal();

            LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(obj.RHS);
            Usol = LHS\RHS;
            obj.Lchol=ichol(LHS);
            obj.L = tril(LHS);
            obj.U = triu(LHS);
            obj.D = diag(diag(LHS));

            obj.computeEigenModes();
            obj.computeModalStiffnessMatrix();
%             obj.createModelPreconditioning();
%             u = obj.solver2(LHS,RHS,refLHS);

            Mid = @(r) r;
            MidOrth = @(r,A,z) z+0.3*(r-A(z));
            LHSf = @(x) LHS*x;
            RHSf = RHS;
          
            %  LHSf = @(x) P*LHS*x;            
            %  RHSf = P*RHS;
            tol = 1e-8;
            P = @(r) Mid(r); %obj.multiplePrec(r,Mid,Mid,LHSf);
            [uCG,residualCG,errCG,errAnormCG] = obj.preconditionedConjugateGradient(LHSf,RHSf,Usol,P,tol);
%             [uCG,residualCG,errCG,errAnormCG] = obj.preconditionedRichardson(LHSf,RHSf,Usol,Mid);

            
            Meifem = @(r) obj.solveEIFEM(r);
            MeifemCG = @(r) obj.solveEIFEMCG(r);
            MeifemContinuous = @(r) obj.solveMEIFEMcontinuous(r);
            Milu = @(r) obj.applyILU(r);
            MiluCG = @(r,A) obj.ILUCG(r,A);
            Mmodal = @(r) obj.solveModalProblem(r);    
            Milu_m = @(r) Milu(Mmodal(r));
            MgaussSeidel = @(r) obj.applyGaussSeidel(r);
            MfixOrth   = @(r,A,z) obj.fixPointOrthogonal(r,A,z);
            
            M = Meifem;%Milu_m;%Meifem; %Milu %Pm
            M2 = MiluCG;
%             [uPCG,residualPCG,errPCG,errAnormPCG] = obj.solverTestEifem(LHSf,RHSf,Usol,M);
            tol = 1e-8;
            P = @(r) obj.multiplePrec(r,M,M2,LHSf);            
            [uPCG,residualPCG,errPCG,errAnormPCG] = obj.preconditionedConjugateGradient(LHSf,RHSf,Usol,P,tol);
%             [uPCG,residualPCG,errPCG,errAnormPCG] = obj.preconditionedRichardson2(LHSf,RHSf,Usol,M,M2);


            figure
            plot(residualPCG,'linewidth',2)
            hold on
            plot(residualCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend({'CG + EIFEM+ ILU(CG-90%-L2)','CG'},'FontSize',12)
            xlabel('Iteration')
            ylabel('Residual')

             figure
            plot(errPCG,'linewidth',2)
            hold on
            plot(errCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            xlabel('Iteration')
            ylabel('||error||_{L2}')

            figure
            plot(errAnormPCG,'linewidth',2)
            hold on
            plot(errAnormCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            xlabel('Iteration')
            ylabel('Energy norm')



            u = obj.solver2(LHS,RHS,refLHS);

% %             obj.obtainCornerNodes();
% %             fineMesh = MeshFromRVE
%             obj.createSubDomainMeshes();
%             obj.createInterfaceSubDomainMeshes();
%             obj.createDomainMesh();

            
%             s.referenceMesh = obj.referenceMesh;
%             mC = MeshCreatorFromSubmeshes();
%             obj.meshDomain = mC.mesh;

%             preconditioner = obj.createPreconditioner(mC.submeshes);

        end

        
        function solveDomainProblem(obj)
            s.mesh     = obj.meshDomain;
            s.bc       = obj.boundaryConditions;
            s.material = obj.material;
            s.type     = 'ELASTIC';
            s.scale    = 'MACRO';
            s.dim      = '2D';
            s.solverTyp = 'ITERATIVE';
            s.iterativeSolverTyp = 'PCG';
            s.preconditionerType = 'EIFEM';
            s.tol = 1e-6;
            
            fem        = FEM.create(s);
            fem.solve();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [15 1]; %nx ny
            obj.scale        = 'MACRO';
            obj.ndimf        = 2;
            obj.functionType = 'P1';
            obj.solverCase   = 'REDUCED';
%             obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
          obj.EIFEMfilename = 'DEF_Q4porL_1.mat';
            obj.weight       = 0.5;
        end

        function createReferenceMesh(obj)
%             filename   = 'lattice_ex1';
%             a.fileName = obj.EIFEMfilename;
%             femD       = FemDataContainer(a);
%             mS         = femD.mesh;
%             bS         = mS.createBoundaryMesh();
%              % Generate coordinates
%             x1 = linspace(0,1,5);
%             x2 = linspace(0,1,5);
%             % Create the grid
%             [xv,yv] = meshgrid(x1,x2);
%             % Triangulate the mesh to obtain coordinates and connectivities
%             [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');
% 
%             s.coord    = coord(:,1:2);
%             s.connec   = F;
%             mS         = Mesh.create(s);
%             bS         = mS.createBoundaryMesh();
            load(obj.EIFEMfilename);
            s.coord    = EIFEoper.MESH.COOR;
            s.coord(s.coord==min(s.coord)) = round(s.coord(s.coord==min(s.coord)));
            s.coord(s.coord==max(s.coord)) = round(s.coord(s.coord==max(s.coord)));
            s.connec   = EIFEoper.MESH.CN;
            mS         = Mesh.create(s);
            bS         = mS.createBoundaryMesh();
            
            obj.meshReference = mS;
            obj.interfaceMeshReference = bS;
            obj.ninterfaces = length(bS);
        end

         function mCoarse = createCoarseMesh(obj)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh();
            mRVECoarse      = MeshCreatorFromRVE(s);
            [mCoarse,meshSubDomainCoarse,interfaceConnecCoarse] = mRVECoarse.create();
        end



        function cMesh = createReferenceCoarseMesh(obj)            
            xmax = max(obj.meshReference.coord(:,1));
            xmin = min(obj.meshReference.coord(:,1));
            ymax = max(obj.meshReference.coord(:,2));
            ymin = min(obj.meshReference.coord(:,2));
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
            connec = [1 2 3 4];
            connec = [2 3 4 1];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end


        function createDomainMaterial(obj)
%             ngaus = 1;        
            m = obj.meshDomain;                      
            obj.material = obj.createMaterial(m);
        end
        
        function L = computeReferenceMeshLength(obj)
            coord = obj.meshReference.coord;
            Lx = max(coord(:,1));
            Ly = max(coord(:,2));
            L = [Lx Ly];
        end     

        function [Dir,PL] = createRawBoundaryConditions(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isBottom  = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
            isTop  = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir{1}.domain    = @(coor) isLeft(coor) | isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

%             Dir{2}.domain    = @(coor) isRight(coor) ;
%             Dir{2}.direction = [2];
%             Dir{2}.value     = -0.1;

            PL.domain    = @(coor) isTop(coor);
            PL.direction = 2;
            PL.value     = -0.1;
        end        

         function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichletFun = [];
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.meshDomain, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end

            pointload = PointLoad(mesh,PL);
             % need this because force applied in the face not in a point
            pointload.values        = pointload.values/size(pointload.dofs,1);
            fvalues                 = zeros(mesh.nnodes*obj.ndimf,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues                 = reshape(fvalues,obj.ndimf,[])';
            pointload.fun.fValues   = fvalues;

            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
        end

%         function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)
%             dirichlet = obj.createBondaryCondition(boundaryMesh,dirchletBc);
%             pointload = obj.createBondaryCondition(boundaryMesh,newmanBc);
%         end
% 
%         function cond = createBondaryCondition(obj,bM,condition)
%             nbound = length(condition.boundaryId);
%             cond = zeros(1,3);
%             for ibound=1:nbound
%                 ncond  = length(condition.dof(nbound,:));
%                 nodeId= unique(bM{condition.boundaryId(ibound)}.globalConnec);
%                 nbd   = length(nodeId);
%                 condition.value{ibound} = condition.value{ibound}/nbd;
%                 for icond=1:ncond
%                     bdcond= [nodeId, repmat(condition.dof(icond),[nbd,1]), repmat(condition.value(icond),[nbd,1])];
%                     cond=[cond;bdcond];
%                 end
%             end
%             cond = cond(2:end,:);
%         end

%         function material = createMaterial(obj,mesh)
%             I = ones(mesh.nelem,obj.quad.ngaus);
%             s.ptype = 'ELASTIC';
%             s.pdim  = '2D';
%             s.nelem = mesh.nelem;
%             s.mesh  = mesh;
%             s.kappa = .9107*I;
%             s.mu    = .3446*I;
%             mat = Material.create(s);
%             mat.compute(s);
%             material = mat;
%         end
        
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

        function dim = getFunDims(obj)
            d.ndimf  = obj.displacementFun.ndimf;
            d.nnodes = size(obj.displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.meshDomain.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

%         function computeStiffnessMatrix(obj)
%             s.type     = 'ElasticStiffnessMatrix';
%             s.mesh     = obj.meshDomain;
%             s.fun      = obj.displacementFun;
%             s.material = obj.material;
%             lhs = LHSintegrator.create(s);
%             obj.LHS = lhs.compute();
%         end

        function LHS = computeStiffnessMatrix(obj,mesh,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.test     = LagrangianFunction.create(s.mesh,obj.ndimf, obj.functionType);
            s.trial    = dispFun;
            s.material = mat;
            s.quadratureOrder = 2;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

%         function computeForces(obj)
%             s.type = 'Elastic';
%             s.scale    = obj.scale;
%             s.dim      = obj.getFunDims();
%             s.BC       = obj.boundaryConditions;
%             s.mesh     = obj.meshDomain;
%             s.material = obj.material;
% %             s.globalConnec = obj.displacementField.connec;
%             s.globalConnec = obj.meshDomain.connec;
%             RHSint = RHSintegrator.create(s);
%             rhs = RHSint.compute();
%             R = RHSint.computeReactions(obj.LHS);
% %             obj.variables.fext = rhs + R;
%             obj.RHS = rhs;
%         end

        function forces = computeForces(obj,stiffness)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
%             s.dim       = obj.getFunDims();
            s.dim.ndofs = obj.displacementFun.nDofs;
            s.BC        = obj.boundaryConditions;
            s.mesh      = obj.meshDomain;
            s.material  = obj.material;
%             s.globalConnec = mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            % Perhaps move it inside RHSint?
            if strcmp(obj.solverCase,'REDUCED')
                R = RHSint.computeReactions(stiffness);
                forces = rhs+R;
            else
                forces = rhs;
            end
            F = obj.assembleVec();
%             forces = forces+F;
        end

        function rhs = assembleVec(obj)
            f = -0.01*ones(obj.displacementFun.nDofsElem,1);
            f(1:2:end) = 0;
            F = repmat(f,[1,1,obj.meshDomain.nelem]);
            s.fun  = obj.displacementFun; % !!!
            assembler = AssemblerFun(s);
            rhs = assembler.assembleV(F, s.fun);
        end

         function  localGlobalDofConnec = createlocalGlobalDofConnec(obj)
            ndimf = obj.displacementFun.ndimf;
            ndom  = obj.nSubdomains(1)*obj.nSubdomains(2);
            for dom = 1:ndom
                row = ceil(dom/obj.nSubdomains(1));
                col = dom-(row-1)*obj.nSubdomains(1);
                localGlobalDofConnecDom = zeros(1,2);
                nodeG = obj.locGlobConnec(:,1,dom);
                nodeL = obj.locGlobConnec(:,2,dom);
                for iunkn = 1:ndimf
                    dofConec = [ndimf*(nodeG - 1) + iunkn ,  ndimf*(nodeL - 1) + iunkn] ;
                    localGlobalDofConnecDom = [localGlobalDofConnecDom;dofConec ];
                end
                localGlobalDofConnec{row,col} = localGlobalDofConnecDom(2:end,:);
            end
         end

         function [interfaceDof,interfaceDom] = computeLocalInterfaceDof(obj)
            intConec = reshape(obj.interfaceConnec',2,obj.interfaceMeshReference{1}.mesh.nnodes,[]);
            intConec = permute(intConec,[2 1 3]);
            nint = size(intConec,3);
            globaldof=0;
            ndimf = obj.ndimf;
            ndofs = obj.meshReference.nnodes*ndimf;
            
            for iint=1:nint
                ndom = size(intConec,2); %length(intConec(1,:,iint));
                for idom = 1:ndom
                    dofaux=0;
                    nodesI = intConec(:,idom,iint);
                    dom = ceil(intConec(1,idom,iint)/obj.meshReference.nnodes);
                    globaldof = (dom-1)*ndofs;
                    for iunkn=1:ndimf
                        DOF = ndimf*(nodesI - 1) + iunkn;
                        DOF = DOF-globaldof;
                        dofaux= [dofaux; DOF];
                    end
                    interfaceDof(:,idom,iint) = dofaux(2:end);
                    interfaceDom(iint,idom) = dom;
                    %                     globaldof = globaldof + (iint*(idom-1)+iint)*dim.ndofs;
                end
                %                 interfaceDof(:,iint) = dofaux(2:end);
                %                 globaldof = globaldof + dim.ndofs;
            end
        end
         
        function Gvec = local2global(obj,Lvec)
%             ndimf  = obj.displacementFun.ndimf;
            Gvec   = zeros(obj.displacementFun.nDofs,obj.nSubdomains(1)*obj.nSubdomains(2));
%             Gvec(locGlobConnec(:,1)) = Lvec(locGlobConnec(:,2));
            ind    = 1;
            for jdom = 1: obj.nSubdomains(2)
                for idom = 1: obj.nSubdomains(1)
                    locGlobConnec = obj.localGlobalDofConnec{jdom,idom};
                    Gvec(locGlobConnec(:,1),ind) = Lvec(locGlobConnec(:,2),ind);
                    ind=ind+1;
                end
            end
        end

        function Lvec = global2local(obj,Gvec)
            ndimf  = obj.displacementFun.ndimf;
            Lvec   = zeros(obj.meshReference.nnodes*ndimf,obj.nSubdomains(1)*obj.nSubdomains(2));
            ind    = 1;
            for jdom = 1: obj.nSubdomains(2)
                for idom = 1: obj.nSubdomains(1)
                    locGlobConnec = obj.localGlobalDofConnec{jdom,idom};
                    Lvec(locGlobConnec(:,2),ind) = Gvec(locGlobConnec(:,1));
                    ind=ind+1;
                end
            end
        end

        function eifem = createEIFEM(obj,meshDomainCoarse,Dir,Kfine)
            RVE         = TrainedRVE(obj.EIFEMfilename);
            s.RVE      = RVE;
            s.mesh     = meshDomainCoarse;
            s.DirCond  = Dir;
            s.Kfine    = Kfine;
            eifem      = EIFEM(s);
        end

        function uInt = computeInterfaceDisp(obj,u)
            nint = size(obj.interfaceDof,3);
            uInt = zeros(size(obj.interfaceDof,1),nint);
            w = [obj.weight,1- obj.weight];
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    unodal = u(:,dom);
                    dof = obj.interfaceDof(:,idom,iint);
                    uInt(:,iint) = uInt(:,iint) + w(idom)*unodal(dof);
                end
            end
        end

        function u = smoothDisplacement(obj,u,uInterface)
            uG = obj.local2global(u);
            uG = sum(uG,2);
            u  = obj.updateInterfaceValues(uG,uInterface); 
        end

        function u = updateInterfaceValues(obj,u,uInterface)
            nint = size(obj.interfaceDof,3);
            for iint = 1:nint
                dom = obj.interfaceDom(iint,1);
                dof = obj.interfaceDof(:,1,iint);
                row = ceil(dom/obj.nSubdomains(1));
                col = dom-(row-1)*obj.nSubdomains(1);
                locGlobConnec = obj.localGlobalDofConnec{row,col};
                [~,ind]    = ismember(dof,locGlobConnec(:,2));
                dofGlob    = locGlobConnec(ind,1);
                u(dofGlob) = uInterface(:,iint);
            end
        end

        function R = scaleInterfaceValues(obj,R)
            nint = size(obj.interfaceDof,3);
            uInt = zeros(size(obj.interfaceDof,1),nint);
            w = [obj.weight,1- obj.weight];
            for iint = 1:nint
                ndom = size(obj.interfaceDof(:,:,iint),2);
                for idom = 1:ndom
                    dom = obj.interfaceDom(iint,idom);
                    Rnodal = R(:,dom);
                    dof = obj.interfaceDof(:,idom,iint);
                    R(dof,dom) = w(idom)* R(dof,dom);
                end
            end
        end

        function Rsbd = computeSubdomainResidual(obj,R)
            RG    = obj.bcApplier.reducedToFullVectorDirichlet(R);
%             obj.plotSolution(RG,obj.meshDomain,0,1,iter,1)
            RGsbd = obj.global2local(RG);
            Rsbd = obj.scaleInterfaceValues(RGsbd);
        end

        function fc = computeContinousField(obj,f)
            fInt  = obj.computeInterfaceDisp(f);
            fc    = obj.smoothDisplacement(f,fInt);
            fc    = obj.bcApplier.fullToReducedVectorDirichlet(fc);
        end

        function z = applyILU(obj,r)
            Lchol=obj.Lchol;
            z = Lchol\r;
            z = (Lchol')\z;
        end

%          function z = applyGaussSeidel(obj,r)
%             L=obj.L;
%             z = L\r;
%              z = (L')\z;
%         end

        function z = applyGaussSeidel(obj,r)
            L=obj.L;
            U=obj.U;
            D=obj.D;
            z = U*r;
            z = D\z;
            z = L*z;
%             z = L\r;
%              z = (L')\z;
        end

        function u = solver(obj,LHS,RHS)
            tol=1e-8;
            iter=1;
            uN = zeros(length(RHS),1);
            R  = RHS - LHS*uN;
            Rsbd = obj.computeSubdomainResidual(R,iter);
            uSbd  = obj.EIFEMsolver.apply(Rsbd);
%             for jdom = 1:obj.nSubdomains(2)
%                 for idom =1:obj.nSubdomains(1)
%                     mesh = obj.meshSubDomain{jdom,idom};
%                     ind = obj.nSubdomains(2)*(jdom-1)+idom;
%                     x = uSbd(:,ind);
%                     row = jdom;
%                     col = idom;
%                     obj.plotSolution(x,mesh,row,col,iter,0)
%                 end
%             end
            u = obj.computeContinousField(uSbd);
            uplot = obj.bcApplier.reducedToFullVectorDirichlet(u);
            obj.plotSolution(uplot,obj.meshDomain,0,11,iter,0)
            uN    = u;
            R  = RHS - LHS*uN;
            Rsbd = obj.computeSubdomainResidual(R,iter);
            uDefSbd  = obj.EIFEMsolver.applySubdomainNeumannDeformational(Rsbd);
            uDef = obj.computeContinousField(uDefSbd);
            uRBsbd = obj.EIFEMsolver.applySubdomainNeumannRigidBody(Rsbd);
            uRB = obj.computeContinousField(uRBsbd);
            uplot = obj.bcApplier.reducedToFullVectorDirichlet(uDef);
%             obj.plotSolution(uplot,obj.meshDomain,0,11,iter+1,0)
            uN    = uN+uDef+uRB;
            uplot = obj.bcApplier.reducedToFullVectorDirichlet(uN);
            obj.plotSolution(uplot,obj.meshDomain,0,11,iter+1,0)
            R  = RHS - LHS*uN;
            e(iter) = norm(R);
            theta = 0.5;
            while e(iter)>tol

                RG    = obj.bcApplier.reducedToFullVectorDirichlet(R);
                obj.plotSolution(RG,obj.meshDomain,0,1,iter,1)
                RGsbd = obj.global2local(RG);
                RGsbd = obj.scaleInterfaceValues(RGsbd);
                uSbd  = obj.EIFEMsolver.apply(RGsbd);
                %                 for jdom = 1:obj.nSubdomains(2)
                %                     for idom =1:obj.nSubdomains(1)
                %                         mesh = obj.meshSubDomain{jdom,idom};
                %                         ind = obj.nSubdomains(2)*(jdom-1)+idom;
                %                         x = uSbd(:,ind);
                %                         row = jdom;
                %                         col = idom;
                %                         obj.plotSolution(x,mesh,row,col,iter,0)
                %                     end
                %                 end
                uInt  = obj.computeInterfaceDisp(uSbd);
                u     = obj.smoothDisplacement(uSbd,uInt);
                u     = obj.bcApplier.fullToReducedVectorDirichlet(u);

                uN    = uN + (1-theta)*u;
                uplot = obj.bcApplier.reducedToFullVectorDirichlet(uN);
                obj.plotSolution(uplot,obj.meshDomain,0,1,iter,0)
                R     = RHS - LHS*uN;
                iter = iter+1;
                e(iter) = norm(R);

            end
        end

         function u = solver2(obj,LHS,RHS,refLHS)
            pinvLHSref = pinv(full(refLHS));
            D = diag(diag(LHS));
            L = tril(LHS);
            U = triu(LHS,1);
%             pinvLHSref = inv(refLHS'*refLHS);
            tol=1e-8;
            iter=1;
            uN = zeros(length(RHS),1);
%             R  = RHS - LHS*uN;
%             Rsbd = obj.computeSubdomainResidual(R,iter);
%             uSbd  = obj.EIFEM.apply(Rsbd);
%             u = obj.computeContinousField(uSbd);
% 
%             uplot = obj.bcApplier.reducedToFullVectorDirichlet(u);
%             obj.plotSolution(uplot,obj.meshDomain,0,11,iter,0)
%             uN    = u;
%             R  = RHS - LHS*uN;
%             Rsbd = obj.computeSubdomainResidual(R,iter);
%             uDefSbd  = obj.EIFEM.applySubdomainNeumannDeformational(Rsbd);
%             uDef = obj.computeContinousField(uDefSbd);
%             uRBsbd = obj.EIFEM.applySubdomainNeumannRigidBody(Rsbd);
%             uRB = obj.computeContinousField(uRBsbd);
%             uplot = obj.bcApplier.reducedToFullVectorDirichlet(uDef);
% %             obj.plotSolution(uplot,obj.meshDomain,0,11,iter+1,0)
%             uN    = uN+uDef+uRB;
%             uplot = obj.bcApplier.reducedToFullVectorDirichlet(uN);
%             obj.plotSolution(uplot,obj.meshDomain,0,11,iter+1,0)
            R  = RHS - LHS*uN;
            e(iter) = norm(R);
            theta = 0.8;
            while e(iter)>tol
                RG    = obj.bcApplier.reducedToFullVectorDirichlet(R);
%                 obj.plotSolution(RG,obj.meshDomain,0,1,iter,1)
%                 Rsbd = obj.computeSubdomainResidual(R,iter);
%                 upinv = pinvLHSref*Rsbd;
%                 upinv = obj.computeContinousField(upinv);
%                 uN    = uN+upinv;
                uplot = obj.bcApplier.reducedToFullVectorDirichlet(uN);
%                  obj.plotSolution(uplot,obj.meshDomain,0,1,iter,0)
%                 R     = RHS - LHS*uN;
                Rsbd = obj.computeSubdomainResidual(R,iter);
                uSbd  = obj.EIFEMsolver.apply(Rsbd);
                u = obj.computeContinousField(uSbd);
                
%                 Rsbd = obj.computeSubdomainResidual(R,iter);
%                 uDefSbd  = obj.EIFEM.applySubdomainNeumannDeformational(Rsbd);
%                 uDef = obj.computeContinousField(uDefSbd);
%                 uRBsbd = obj.EIFEM.applySubdomainNeumannRigidBody(Rsbd);
%                 uRB = obj.computeContinousField(uRBsbd);

%                 uN    = uN + (1-theta)*(u+uDef+uRB);
                uN1    = uN + (u);
                R     = RHS - U*uN1;
%                 z = obj.applyILU(R);
%                 uN1    = uN1 + z;
%                 R     = RHS - U*uN1;
                uN    = L\(R);
%                 uN    = theta*uN + (1-theta)*L\(R);
%                 uplot = obj.bcApplier.reducedToFullVectorDirichlet(uN);
%                 obj.plotSolution(uplot,obj.meshDomain,0,1,iter,0)
                R  = RHS - LHS*uN;
                iter = iter+1;
                e(iter) = norm(R);

            end
         end

         function [x,residual,err,errAnorm] = preconditionedRichardson(obj,A,B,xsol,P)
            tol = 1e-8;
            iter = 0;
            n = length(B);
            x = zeros(n,1);
            r = B - A(x);      
            z = P(r);
%             x = x + z;
%             r = B - A(x);  
            theta = 0.999;
            tau =0.1;
            while norm(r) > tol
                x = x + tau * z;
                r = B - A(x); 
                z = P(r);
                test(iter+1)=norm(z);
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A(x-xsol);                
            end
         end   

         function [x,residual,err,errAnorm] = preconditionedRichardson2(obj,A,B,xsol,P)
            tol = 1e-8;
            iter = 0;
            n = length(B);
            x = zeros(n,1);
            r = B - A(x);      
            z = P(r);
            theta = 0.1;
%             tau =(r'*z)/(z'*A(z));
            tau=0.2;
            while norm(r) > tol
                x = x + tau * z;
                r = B - A(x); 
                z = P(r); 
                test(iter+1)=norm(z);
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A(x-xsol);                
            end
         end   
         
         function [x,residual,err,errAnorm] = preconditionedConjugateGradient(obj,A,B,xsol,P,tol)
            iter = 0;
            n = length(B);
            x = zeros(n,1);
            r = B - A(x);               
            z = P(r);
            p = z;
            rzold = r' * z; 
            while norm(r) > tol
                Ap = A(p);
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                z = P(r);
                rznew = r' * z;
                p = z + (rznew / rzold) * p;
                rzold = rznew;
                iter = iter + 1;
                residual(iter) = norm(r); 
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A(x-xsol);                
            end
         end



         function [x,residual,err,errAnorm] = solverTestEifem(obj,A,B,xsol,P)
            tol = 1e-8;
            iter = 0;
            n = length(B);
            x = zeros(n,1);
%             load("xEifem.mat");
            r = B - A(x);    
            RGsbd = obj.computeSubdomainResidual(r);
            uSbd =  obj.EIFEMsolver.apply(RGsbd);
            z = P(r);
            p = z;
            rzold = r' * z;           
            while norm(r) > tol
                Ap = A(p);
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                z = P(r);
%                 z=r;
                rznew = r' * z;
                p = z + (rznew / rzold) * p;
                rzold = rznew;
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A(x-xsol);                
            end
         end

         function computeEigenModes(obj)
            LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
            nbasis = 8;
            [phi,D]=eigs(LHS,nbasis,'smallestabs');
            obj.eigenModes = phi;
         end

         function computeModalStiffnessMatrix(obj)
            LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);             
            phi = obj.eigenModes;
            lhs = phi'*LHS*phi;
            obj.Kmodal = lhs;
         end

         function createModelPreconditioning(obj)
            phi=obj.eigenModes;
            LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);             
            I = speye(size(LHS,1));
            M = I-phi*((phi'*LHS*phi)\(phi'));
            obj.MmodalPrecond = M;
         end

        

        function z = solveModalProblem(obj,r)
             lhs=obj.Kmodal;
             phi=obj.eigenModes;
             LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);  
             r1=phi'*r;
             zP=lhs\r1;
             z =phi*zP;
%               z = (r-LHS*z);
             
            %M = obj.MmodalPrecond;
            %z = M*r;
        end  
     
    
        function z = solveEIFEM(obj,r)
            RGsbd = obj.computeSubdomainResidual(r);
            uSbd =  obj.EIFEMsolver.apply(RGsbd);
            z = obj.computeContinousField(uSbd);
            LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);  
%             z = r - LHS*z;  
%             alpha =    0.659336578745820;% rN'*rN/(rN'*LHS*rN);
%             z= r + z;
        end

        function z = solveEIFEMCG(obj,r)
%             xsol = obj.bcApplier.reducedToFullVectorDirichlet(xsol);
            tol = 1e-8;
            iter = 0;
            RGsbd = obj.computeSubdomainResidual(r);
            nsbd = size(RGsbd,2);
            r = reshape(RGsbd,[],1);
            n = length(r);
            x = zeros(n,1);
%             load("xEifem.mat");             
            uSbd =  obj.EIFEMsolver.apply(RGsbd);
%             z = r- reshape(obj.refLHS*uSbd,[],1);
            z = reshape(uSbd,[],1);
%             Ax = obj.refLHS*uSbd;
%             Ax = reshape(Ax,[],1);
%             r = r - Ax;  
%             z = obj.solveEIFEM(r);
            p = z;
            rzold = r' * z;   
            rk=r+1;
            while norm(r-rk) > tol
                rk=r;
                pMat = reshape(p,[],nsbd);
                Ap = reshape(obj.refLHS*pMat,[],1);
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                rsbd = reshape(r,[],nsbd);
                uSbd =  obj.EIFEMsolver.apply(rsbd);
%                 z = r- reshape(obj.refLHS*uSbd,[],1);
                z = reshape(uSbd,[],1);
%                 z=r;
                rznew = r' * z;
                p = z + (rznew / rzold) * p;
                rzold = rznew;
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
%                 err(iter)=norm(x-xsol);
%                 errAnorm(iter)=((x-xsol)')*A(x-xsol);                
            end
%              z = r- reshape(obj.refLHS*uSbd,[],1);
%              z=r;
%              p = z;
%             rzold = r' * z;   
%             rk=r+1;
%             while norm(r-rk) > tol
%                 rk=r;
%                 pMat = reshape(p,[],nsbd);
%                 Ap = reshape(obj.refLHS*pMat,[],1);
%                 alpha = rzold / (p' * Ap);
%                 x = x + alpha * p;
%                 r = r - alpha * Ap;
%                 rsbd = reshape(r,[],nsbd);
%                 uSbd =  obj.EIFEMsolver.apply(rsbd);
%                 z = r- reshape(obj.refLHS*uSbd,[],1);
% %                 z=r;
% %                 z = reshape(uSbd,[],1);
% %                 z=r;
%                 rznew = r' * z;
%                 p = z + (rznew / rzold) * p;
%                 rzold = rznew;
%                 iter = iter + 1;
%                 residual(iter) = norm(r); %Ax - b
% %                 err(iter)=norm(x-xsol);
% %                 errAnorm(iter)=((x-xsol)')*A(x-xsol);                
%             end
%             rsbd = reshape(r,[],nsbd);
%             uSbd =  obj.EIFEMsolver.apply(rsbd);
%             z =  reshape(uSbd,[],1);
%             z= z+x;
%             z = reshape(z,[],nsbd);
%             z = obj.computeContinousField(z);  
            z = reshape(x,[],nsbd);
            z = obj.computeContinousField(z);  
        end

        function Ug = computeAssembledProjectionMatrix(obj)
            Ueifem = obj.EIFEMsolver.getProjectionMatrix();
            sizeU = size(Ueifem);
            Ueifem = repmat(Ueifem,1,obj.nSubdomains(1)*obj.nSubdomains(2));
            Ueifem = obj.scaleProjectionMatrix(Ueifem);
            Ug = obj.assembleProjectionMatrix(Ueifem);
        end

        function U = scaleProjectionMatrix(obj,U)
            nDofcoarse = size(U,2)/obj.nSubdomains(1)*obj.nSubdomains(2);
            for i=1: nDofcoarse
                Usbd = U(:,i:nDofcoarse:end);
                Usbd = obj.scaleInterfaceValues(Usbd);
                U(:,i:nDofcoarse:end)= Usbd;
            end
        end
    

        function Ug = assembleProjectionMatrix(obj,U)
            nsbd = obj.nSubdomains(1)*obj.nSubdomains(2);
            nDofcoarse = size(U,2)/nsbd;
            Ug = zeros(obj.displacementFun.nDofs,nDofcoarse*nsbd);
            for i=1: nDofcoarse
                Usbd = U(:,i:nDofcoarse:end);
                Uaux = obj.local2global(Usbd);
                Ug(:,i:nDofcoarse:end)= Uaux;
            end
%             nsbd = obj.nSubdomains(1)*obj.nSubdomains(2);
%             nDofFine = size(U,1);
%             nDofcoarse = size(U,2)/nsbd;
%             Ug = zeros(nDofFine*nsbd,nDofcoarse*nsbd);
%             for i = 1:nsbd
%                 col = nDofcoarse*(i-1)+1;
%                 Usbd = U(:,col:col+nDofcoarse-1);
%                 row = nDofFine*(i-1)+1;
%                 Ug(row:row+nDofFine-1,col:col+nDofcoarse-1) = Usbd;
%             end
            Ug = sparse(Ug);
        end

        function computeKEIFEMglobal(obj)
            U = obj.computeAssembledProjectionMatrix();
%             U = U(:,9:end);
            obj.EIFEMprojection = obj.computeReducedProjection(U);
             LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);  
            obj.KeifemContinuous = obj.EIFEMprojection'*LHS*obj.EIFEMprojection;
            
        end

        function Ured = computeReducedProjection(obj,U)
            ncol = size(U,2);
            for i=1:ncol
                Ured(:,i) = obj.bcApplier.fullToReducedVectorDirichlet(U(:,i));  
            end
        end

        function z = solveMEIFEMcontinuous(obj,r)
             lhs=obj.KeifemContinuous;
             phi=obj.EIFEMprojection;
%              phi = obj.bcApplier.fullToReducedMatrixDirichlet(phi);  
             r1=phi'*r;
             zP=lhs\r1;
             z =phi*zP;
%               z = (r-LHS*z);
             
            %M = obj.MmodalPrecond;
            %z = M*r;
        end  

         function z = multiplePrec(obj,r,P,P2,A)
             z = P(r);
             r=r-A(z);                
             x = P2(r,A);
            z=z+x;                         
         end        

         function x = ILUCG(obj,r,A)
            P = @(r) obj.applyILU(r);
            xsol = zeros(size(r));
            factor = 0.1;
            tol = factor*norm(r);
            [x,residual,err,errAnorm] = obj.preconditionedConjugateGradient(A,r,xsol,P,tol);
         end

        %    P = @(r) obj.multiplePrec(r,Mid,Mid,LHSf);


         function z = fixPointOrthogonal(obj,r,A,zP1)
            r=r-A(zP1);   
            factor = 0.05;
            iter = 0;
            rplot = obj.bcApplier.reducedToFullVectorDirichlet(r);
%             obj.plotSolution(rplot,obj.meshDomain,11,11,iter,1)
            n = length(r);
            x = zeros(n,1);         
            norm0 = norm(r);
            norm0= 1/norm0;
            tau=0.3;
            
            while norm0*norm(r)>factor  
                alpha = r'*r/(r'*A(r));
                x = x + alpha * r;
                r = r - alpha*A(r);
%                 r= r-tau * r;
                iter = iter + 1;     
                rplot = obj.bcApplier.reducedToFullVectorDirichlet(r);
%                 obj.plotSolution(rplot,obj.meshDomain,11,11,iter,1)
            end
            z=zP1+x;
        end
       

        function plotSolution(obj,x,mesh,row,col,iter,flag)
            if nargin <7
                 flag =0;
            end
            %             xFull = bc.reducedToFullVector(x);
            if size(x,2)==1
                 s.fValues = reshape(x,2,[])';
            else
                 s.fValues = x; 
            end
%            
            
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 2;
            s.order = obj.functionType;
            xF = LagrangianFunction(s);
%             xF.plot();
            if flag == 0
                xF.print(['domain',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 1
                xF.print(['DomainResidual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 2
                xF.print(['Residual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 3
                xF.print(['domainFine',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 4
                xF.print(['domainNeuman',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            end
            fclose('all');
        end

    end
end
