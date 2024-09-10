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

        displacementFun
        LHS
        RHS
        scale
        ndimf    
        functionType
        EIFEM
    end

    methods (Access = public)

        function obj = EIFEMtesting()
            close all
            obj.init();
            obj.createReferenceMesh();
           
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.meshReference;         
            m = MeshCreatorFromRVE(s);
            [obj.meshDomain,obj.meshSubDomain,obj.interfaceConnec,~,obj.locGlobConnec] = m.create();

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
           
            cMesh           = createCoarseMesh(obj);
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = cMesh;         
            mRVECoarse      = MeshCreatorFromRVE(s);
            [meshDomainCoarse,meshSubDomainCoarse,interfaceConnecCoarse] = mRVECoarse.create();

            refDisp     = LagrangianFunction.create(obj.meshReference, obj.ndimf,obj.functionType);
            refMat      = obj.createMaterial(obj.meshReference);
            refLHS      = obj.computeStiffnessMatrix(obj.meshReference,refDisp,refMat);

            obj.EIFEM = obj.createEIFEM(meshDomainCoarse,Dir,refLHS);

            LHS = obj.bcApplier.fullToReducedMatrixDirichlet(obj.LHS);
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(obj.RHS);
            Usol = LHS\RHS;
            obj.Lchol=ichol(LHS);
            obj.L = tril(LHS);
            obj.U = triu(LHS);
            obj.D = diag(diag(LHS));

%             u = obj.solver2(LHS,RHS,refLHS);

            [uCG,residualCG,errCG,errAnormCG]  = obj.conjugateGradient(LHS,RHS,Usol);
            [uPCG,residualPCG,errPCG,errAnormPCG] = obj.preconditionedConjugateGradient(LHS,RHS,Usol);

            figure
            plot(residualPCG,'linewidth',2)
            hold on
            plot(residualCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend({'CG + EIFEM','CG'},'FontSize',12)
            xlabel('Iteration')
            ylabel('Residual')

             figure
            plot(errPCG,'linewidth',2)
            hold on
            plot(errCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM','CG')
            xlabel('Iteration')
            ylabel('||error||_{L2}')

            figure
            plot(errAnormPCG,'linewidth',2)
            hold on
            plot(errAnormCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM','CG')
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
            obj.nSubdomains  = [20 1]; %nx ny
            obj.scale        = 'MACRO';
            obj.ndimf        = 2;
            obj.functionType = 'P1';
            obj.solverCase   = 'REDUCED';
%             obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
          obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_8_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_1.mat';
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
            s.connec   = EIFEoper.MESH.CN;
            mS         = Mesh.create(s);
            bS         = mS.createBoundaryMesh();
            
            obj.meshReference = mS;
            obj.interfaceMeshReference = bS;
            obj.ninterfaces = length(bS);
        end

        function cMesh = createCoarseMesh(obj)
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
            eifem       = EIFEM(s);
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

        function Rsbd = computeSubdomainResidual(obj,R,iter)
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
%              z = (L')\z;
        end

        function u = solver(obj,LHS,RHS)
            tol=1e-8;
            iter=1;
            uN = zeros(length(RHS),1);
            R  = RHS - LHS*uN;
            Rsbd = obj.computeSubdomainResidual(R,iter);
            uSbd  = obj.EIFEM.apply(Rsbd);
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
            uDefSbd  = obj.EIFEM.applySubdomainNeumannDeformational(Rsbd);
            uDef = obj.computeContinousField(uDefSbd);
            uRBsbd = obj.EIFEM.applySubdomainNeumannRigidBody(Rsbd);
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
                uSbd  = obj.EIFEM.apply(RGsbd);
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
                uSbd  = obj.EIFEM.apply(Rsbd);
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


        function [x,residual,err,errAnorm] = preconditionedConjugateGradient(obj,A,B,xsol)
            tol = 1e-8;
            iter = 0;
            n = length(B);
            x = zeros(n,1);
            r = B - A * x;
%             r = obj.applyILU(r);
%             r = obj.applyGaussSeidel(r);
%             r = B - A * (x+0.2*z1);
            %             z = ModalTesting.applyPreconditioner(M,r);
            RG  = obj.bcApplier.reducedToFullVectorDirichlet(r);
%             obj.plotSolution(RG,obj.meshDomain,0,1,iter,1)
            RGsbd = obj.global2local(RG);
            RGsbd = obj.scaleInterfaceValues(RGsbd);
            uSbd =  obj.EIFEM.apply(RGsbd);
            uInt  = obj.computeInterfaceDisp(uSbd);
            u     = obj.smoothDisplacement(uSbd,uInt);
            z     = obj.bcApplier.fullToReducedVectorDirichlet(u);
%             z=r-z;
%             z = z-A*obj.applyGaussSeidel(z);
            zilu = obj.applyILU(r);
            z=r-z-zilu;
%             zILU = obj.applyILU(z);
%             z=z-zILU;
            %             z=r-z;
            p = z;
            rzold = r' * z;
            

            while norm(r) > tol
                Ap = A * p;
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
%                 uplot = obj.bcApplier.reducedToFullVectorDirichlet(x);
%               obj.plotSolution(uplot,obj.meshDomain,0,1,iter,0)
                r = r - alpha * Ap;
%                 r = obj.applyILU(r);
%                 r = obj.applyGaussSeidel(r);
%                 r = B - A * (x+0.2*z1);
%                 RG  = obj.bcApplier.reducedToFullVectorDirichlet(r);
% %                  obj.plotSolution(RG,obj.meshDomain,0,1,iter,1)
%                 RGsbd = obj.global2local(RG);
%                 RGsbd = obj.scaleInterfaceValues(RGsbd);
                RGsbd = obj.computeSubdomainResidual(r,iter);
                uSbd =  obj.EIFEM.apply(RGsbd);
%                 uInt  = obj.computeInterfaceDisp(uSbd);
%                 u     = obj.smoothDisplacement(uSbd,uInt);
%                 z     = obj.bcApplier.fullToReducedVectorDirichlet(u);
                z = obj.computeContinousField(uSbd);
                test(iter+1) = norm(z);
                test2(iter+1) = max(z(2:2:end));
                z = r-z;
%                 z = z-A*obj.applyGaussSeidel(z);
                zilu = obj.applyILU(r);
                z=r-z-zilu;
%                 zILU = obj.applyILU(z);
%                 z=z-zILU;
                %                 z = obj.applyPreconditioner(r);
                rznew = r' * z;
                %hasNotConverged = sqrt(rsnew) > tol;


                p = z + (rznew / rzold) * p;
                rzold = rznew;

                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A*(x-xsol);
                
            end
        end

       

        function [x,residu,err,errAnorm] = conjugateGradient(obj,LHS,RHS,xsol)
            tol = 1e-8;
            iter = 1;
            x = zeros(size(RHS));
            r = RHS - LHS * x; 
%             RG  = obj.bcApplier.reducedToFullVectorDirichlet(r);
% %             obj.plotSolution(RG,obj.meshDomain,0,1,iter,1)
%             RGsbd = obj.global2local(RG);
%             RGsbd = obj.scaleInterfaceValues(RGsbd);
%             uSbd =  obj.EIFEM.apply(RGsbd);
%             uInt  = obj.computeInterfaceDisp(uSbd);
%             u     = obj.smoothDisplacement(uSbd,uInt);
%             u     = obj.bcApplier.fullToReducedVectorDirichlet(u);
%             uplot = obj.bcApplier.reducedToFullVectorDirichlet(u);
%             obj.plotSolution(uplot,obj.meshDomain,0,1,iter,0)
%             x = u;
%             r = RHS - LHS * x;
            p = r; 
            rsold = r' * r;
            

            hasNotConverged = true;

            while hasNotConverged
                Ap = LHS * p;
                alpha = rsold / (p' * Ap);
                x = x + alpha * p;
                uplot = obj.bcApplier.reducedToFullVectorDirichlet(x);
%                 obj.plotSolution(uplot,obj.meshDomain,0,1,iter,0)
                r = r - alpha * Ap;
                rsnew = r' * r;

                hasNotConverged = norm(LHS*x - RHS) > tol;

                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
                iter = iter + 1;
                 residu(iter) = norm(LHS*x - RHS); %Ax - b
                 err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*LHS*(x-xsol);
%                 res = LHS*x - RHS;
                
                %conjugateGradient_Solver.plotSolution(x,mesh,bc,iter)
                
                %conjugateGradient_Solver.plotRes(res,mesh,bc,iter)
            end
            %save('residuConjugateZeros.mat', 'residu')
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
