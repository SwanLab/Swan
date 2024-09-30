classdef GeneralPreconditioner < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
             nSubdomains

             Lchol

        L
        U
        D

        scale
                ndimf    
                        functionType
        solverCase



    end

    methods (Access = public)

        function obj = GeneralPreconditioner()
            close all
            obj.init();
            obj.contructionObject()
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
        
        function z = applyILU(obj,r)
            Lchol=obj.Lchol;
            z = Lchol\r;
            z = (Lchol')\z;
        end        

        function x = ILUCG(obj,r,A)
            P = @(r) obj.applyILU(r);
            x0 = zeros(size(r));
            factor = 0.99;
            tol = factor*norm(r);
            x = PCG.solve(A,r,x0,P,tol);
        end
        
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
        


    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [5 1]; %nx ny
            obj.scale        = 'MACRO';
            obj.ndimf        = 2;
            obj.functionType = 'P1';
            obj.solverCase   = 'REDUCED';
            %             obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            obj.EIFEMfilename = 'DEF_Q4porL_1.mat';
            %           obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            obj.weight       = 0.5;
        end

        function constractionObject(obj)
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

            refDisp        = LagrangianFunction.create(obj.meshReference, obj.ndimf,obj.functionType);
            refMat          = obj.createMaterial(obj.meshReference);
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
            Dir{1}.domain    = @(coor) isLeft(coor)| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

            %             Dir{2}.domain    = @(coor) isRight(coor) ;
            %             Dir{2}.direction = [2];
            %             Dir{2}.value     = 0;

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
            w = [obj.weight,1-obj.weight];
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



        %          function z = applyGaussSeidel(obj,r)
        %             L=obj.L;
        %             z = L\r;
        %              z = (L')\z;
        %         end




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

        function z = multiplePrec(obj,r,P1,P2,P3,A)
            z1 = P1(r);
            r  = r-A(z1);
            z2 = P2(r);
            r  = r-A(z2);
            z3 = P3(r);
            z  = z1+z2+z3;

        end

        function z = additivePrec(obj,r,P1,P2,A)
            z1 = P1(r);
            z2 = P2(r);
            z  = z1+z2;

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