classdef CoarseTesting_2Params< handle

    properties (Access = public)
        SolExact
        Sol
        residualPCG
        errPCG
        errAnormPCG
        residualCG
        errCG
        errAnormCG
        residualILU
    end

    properties (Access = private)
        nSubdomains
        tFrame
        tCross
        ic
        icr
        lg
        bs
        meshDomain
        referenceMesh
        subdomainMeshes
        cellMesh
        discMesh
        boundaryConditions
        bcApplier
        ddDofManager

        tolSameNode
        data
        params

        xmin 
        xmax 
        ymin 
        ymax  
        designVariable
        unfittedMesh

    end


    methods (Access = public)

        function obj = CoarseTesting_2Params(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.createMesh();
            mR  = obj.referenceMesh;
            bS  = mR.createBoundaryMesh();                         
            [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
            obj.boundaryConditions = bC;
            obj.createBCapplier()
            [LHS,RHS,~] = obj.createElasticProblem();

            %EXACT SOLUTION
            tic
            LHSf   = @(x) LHS*x;
            RHSf   = RHS;
            Usol   = LHS\RHS;
            Ufull  = obj.bcApplier.reducedToFullVectorDirichlet(Usol); 
            t_direct=toc

            % PRECONDITIONERS
            Milu        = obj.createILUpreconditioner(LHS);
            Mid         = @(r)r;
            Mcoarse     = obj.createCoarseNNPreconditioner(mR,dir,obj.ic,obj.lg,bS,obj.icr,obj.discMesh);
            Mmult        = @(r) Preconditioner.multiplePrec(r,LHSf,Milu,Mcoarse,Milu);


            tol = 1e-8;
            x0  = zeros(size(RHSf));

            tic  %SOLVE THE CASE WITH STANDARD CG
            [~,obj.residualCG,errCG, errCG] = PCG.solve(LHSf,RHSf,x0,Mid,tol,Usol,obj.meshDomain,obj.bcApplier);
            t_CG=toc
            tic  % SOLVE THE CASE WITH CG+ ILU
            [~,obj.residualILU,errILU, errAnormILU] = PCG.solve(LHSf,RHSf,x0,Milu,tol,Usol,obj.meshDomain,obj.bcApplier);
            t_ILU=toc
            tic  % SOLVE THE CASE WITH PRECONDITIONING ILU+EIFEM+ILU
            [uPCG,obj.residualPCG,obj.errPCG,obj.errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol,obj.meshDomain,obj.bcApplier);
            t_PCG=toc
            xFull = obj.bcApplier.reducedToFullVectorDirichlet(uPCG);
            

            uDomain = obj.bcApplier.reducedToFullVectorDirichlet(uPCG);
            % uDomain = obj.ddDofManager.global2local(uDomain);
            
            % LAGRANGIAN FUN SOLUTIONS
            s.mesh     = obj.meshDomain;
            s.ndimf    = obj.meshDomain.ndim;
            s.order    = 'P1';
            s.fValues  = reshape(xFull,2,[])';
            obj.Sol    = LagrangianFunction(s); %Preconditioned sol  
            s.fValues = reshape(Ufull,2,[])';
            obj.SolExact=LagrangianFunction(s); %Exact sol

            % obj.print(uPCG,"SolPCG");
            %CoarsePlotSolution(uFun, obj.meshDomain, obj.bcApplier,'TestCoarseAbril', obj.r, obj.centroids);
            %CoarsePlotSolution(RealFun, obj.meshDomain, obj.bcApplier,'TestRealAbril', obj.r, obj.centroids);

        end

        function PlotSolution(obj)
            figure('Position',[65 280 1800 600])
            tiledlayout(1,3)

            nexttile
            plot(obj.residualPCG,'linewidth',2)
            hold on
            plot(obj.residualILU,'linewidth',2)
            set(gca, 'YScale', 'log')
            xlabel('Iteration')
            ylabel('Residual')
            title("Residual")
            legend({'PCG','CG'})

            nexttile
            plot(obj.errPCG,'linewidth',2)
            hold on
            plot(errILU,'linewidth',2)
            set(gca, 'YScale', 'log')
            xlabel('Iteration')
            ylabel('||error||_{L2}')
            title("error")
            legend({'PCG','CG'})

            nexttile
            plot(obj.errAnormPCG,'linewidth',2)
            hold on
            plot(errAnormILU,'linewidth',2)
            set(gca, 'YScale', 'log')
            xlabel('Iteration')
            ylabel('Energy norm')
            title("Err Anorm")
            legend({'PCG','CG'})
        end

        function print(obj,sol,fileName)
                z.mesh      = obj.meshDomain;
                z.order     = 'P1';
                z.fValues   = reshape(sol,z.mesh.ndim,[])';
                uFeFun = LagrangianFunction(z);%
                uMeshFun = obj.unfittedMesh.obtainFunctionAtUnfittedMesh(uFeFun);
                
                fvalues = [uMeshFun.innerMeshFunction.fValues;
                    uMeshFun.innerCutMeshFunction.fValues];

                s.coord = [uMeshFun.innerMeshFunction.mesh.coord;
                    uMeshFun.innerCutMeshFunction.mesh.coord];

                s.connec = [uMeshFun.innerMeshFunction.mesh.connec;
                    uMeshFun.innerCutMeshFunction.mesh.connec  + max(uMeshFun.innerMeshFunction.mesh.connec(:))];
                
                mh = Mesh.create(s);
                ss.mesh = mh;
                ss.fValues = fvalues;
                ss.order = 'P1';
                ss.ndimf = size(fvalues,2)
                u = LagrangianFunction(ss);
                u.print(fileName,'Paraview')
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        function init(obj,cParams)
            p.Training      = cParams.Training;    % 'EIFEM'/'Multiscale'/'EIFisol'
            p.Option        = cParams.Option;      % 'Dataset'/'NN'
            p.nelem         = cParams.nelem;       %  Mesh refining
            obj.params      = p;
            obj.tFrame      = cParams.tFrame;
            obj.tCross      = cParams.tCross;
            obj.nSubdomains = size(obj.tFrame');
            obj.tolSameNode = 1e-11;
        end

        function createMesh(obj)
            mSbd = obj.createSubDomainMeshes();
            [mD,mSb,iC,lG,iCR,discmesh] = obj.createMeshDomainJoiner(mSbd);
            obj.meshDomain = mD;        % mD:conj subdominis --> Tot el domini
            obj.subdomainMeshes = mSb;  %??? % mSb: subdonain Meshes
            obj.ic              = iC;   % interface Connectivities ???
            obj.icr             = iCR;  % info de les coordenades del corresponent subdomini 
            obj.lg              = lG;   % localGlobal 
            obj.discMesh        = discmesh;
        end


        function  mSbd = createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            Lx = 2; Ly = 2;
            cM=cell(nY,nX);
            mSbd=cell(nY,nX);
            for jDom = 1:nY
                for iDom = 1:nX
                    refMesh  =obj.createStructuredMesh();
                    coord0  = refMesh.coord;
                    s.coord(:,1) = coord0(:,1)+Lx*(iDom-1);
                    s.coord(:,2) = coord0(:,2)+Ly*(jDom-1);
                    s.connec = refMesh.connec;
                    mIJ     = Mesh.create(s);
                    %                     plot(mIJ)
                    %                     hold on;
                    mSbd{jDom,iDom} = mIJ;
                    cM{jDom,iDom} = refMesh;  %same but with local coordinates
                end
            end
            obj.referenceMesh = mSbd{1,1};
            obj.cellMesh=cM;   
        end


        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomainJoiner(obj,mSbd)
           s.nsubdomains   = obj.nSubdomains; %nx ny
           s.meshReference = obj.referenceMesh;
           s.tolSameNode   = obj.tolSameNode;
           s.meshSbd       = mSbd;
           m = MeshJoiner(s);
           [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end

        function mS = createStructuredMesh(obj)
            n =obj.params.nelem;
            x1      = linspace(-1,1,n);
            x2      = linspace(-1,1,n);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            % 
            % mesh= QuadMesh(1,1,n,n);
            % s.coord= mesh.coord;
            % s.connec=mesh.connec;

            obj.xmin = min(x1);            
            obj.xmax = max(x1);
            obj.ymin = min(x2);
            obj.ymax = max(x2);
            delta = 1e-8;
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-delta];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[delta,-delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[delta,delta];
            mS = Mesh.create(s);
        end
        
        function mCoarse = createCoarseMesh(obj)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh();
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE2D(s);
            [mCoarse,~,~] = mRVECoarse.create();
        end


        function cMesh = createReferenceCoarseMesh(obj)
            coord(1,1) = obj.xmin;        % Crea els nodes i assigna als DOfs
            coord(1,2) = obj.ymin;        % la coordenada corresponent
            coord(2,1) = obj.xmax;
            coord(2,2) = obj.ymin;
            coord(3,1) = obj.xmax;
            coord(3,2) = obj.ymax;
            coord(4,1) = obj.xmin;
            coord(4,2) = obj.ymax;

            connec = [1 2 3 4];    % crea conectivitats entre els 4 nodes
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);  % crea la mesh de 4 nodes
        end

        function createBCapplier(obj)
            s.mesh                  = obj.meshDomain;
            s.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(s);
        end


        function material = createMaterial(obj)
            obj.createDesignVariable()
            m= obj.createMaterialInterpolator();
            s.type                 = 'DensityBased';
            s.density              = obj.designVariable;
            s.materialInterpolator = m;
            s.dim                  = '2D';
            s.mesh                 = obj.meshDomain;
            material = Material.create(s);
            material.setDesignVariable({obj.designVariable.fun})
            material = material.obtainTensor();
        end


        function createDesignVariable(obj)
            ls= obj.computeLevelSet();
            sUm.backgroundMesh = obj.meshDomain;
            sUm.boundaryMesh   = obj.meshDomain.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(-ls);
            % aaa=uMesh.createInnerMesh();
            Mprint=UnfittedMesh(sUm);
            Mprint.compute(ls);
            obj.unfittedMesh=Mprint;
            funLS        = CharacteristicFunction.create(uMesh);
            s.filterType = 'LUMP';
            s.mesh       = obj.meshDomain;
            s.trial      = LagrangianFunction.create(obj.meshDomain,1,'P1');
            f = Filter.create(s);
            s.fun      = f.compute(funLS,2);
            s.type     = 'Density';
            s.plotting = false;
            dens = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function ls=computeLevelSet(obj)
            [x0,y0] = obj.computeSubdomainCentroid();
            [Nx,Ny] = size(obj.tFrame);
            GeomParams(Nx,Ny) = struct('type',[],'length',[],'tFrame',[],'tCross',[],'xCoorCenter',[],'yCoorCenter',[]);

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    GeomParams(i,j).type        = "CrossedSquare";
                    GeomParams(i,j).length      = 2;
                    GeomParams(i,j).tFrame      = obj.tFrame(i,j);
                    GeomParams(i,j).tCross      = obj.tCross(i,j);
                    GeomParams(i,j).xCoorCenter = x0(i,j);
                    GeomParams(i,j).yCoorCenter = y0(i,j);
                end
            end
            s.type        = 'GivenPattern';
            s.paramsList  = GeomParams;
            g             = GeometricalFunction(s);
            phiFun        = g.computeLevelSetFunction(obj.meshDomain);
            ls            = phiFun.fValues;
         end

         function [x0,y0]= computeSubdomainCentroid(obj)
            % [Nx, Ny] = size(obj.r);
            xMin=min(obj.meshDomain.coord(:,1));
            xMax=max(obj.meshDomain.coord(:,1));
            yMin=min(obj.meshDomain.coord(:,2));
            yMax=max(obj.meshDomain.coord(:,2));

            dx = obj.xmax-obj.xmin;
            dy = obj.ymax-obj.ymin;

            x_center = xMin + dx/2 : dx : xMax - dx/2;
            y_center = yMin + dy/2 : dy : yMax - dy/2;

            [x0, y0] = meshgrid(x_center, y_center);
        end

        function m=createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.meshDomain.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            E1 = 1;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA           = matA;
            s.matB           = matB;

            m = MaterialInterpolator.create(s);
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            minx = min(obj.meshDomain.coord(:,1));
            maxx = max(obj.meshDomain.coord(:,1));
            miny = min(obj.meshDomain.coord(:,2));
            maxy = max(obj.meshDomain.coord(:,2));
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            isBottom = @(coor) (abs(coor(:,2) - miny)   < tolBound);
            isTop    = @(coor) (abs(coor(:,2) - maxy)   < tolBound);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

            PL.domain    = @(coor) isRight(coor);
            PL.direction = 2;
            PL.value     = -1;       %Set displacement intensity 
        end 

        function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichletFun = [];
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.meshDomain, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end
            pointload = TractionLoad(mesh,PL,'DIRAC');
            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
        end


        function [LHSr,RHSr,lhs] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial();
            [lhs,LHSr] = obj.computeStiffnessMatrix(u,obj.meshDomain,material);
            RHSr       = obj.computeForces(lhs,u);
        end

        function [LHS,LHSr] = computeStiffnessMatrix(obj,dispFun,mesh,C)
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            LHS= IntegrateLHS(f,dispFun,dispFun,mesh,'Domain',2);
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

        function Mcoarse = createCoarseNNPreconditioner(obj,mR,dir,iC,lG,bS,iCR,dMesh)
            p=obj.params;
            switch p.Option
                case 'Dataset'
                    nameFile=obj.computeNameFile();
                    obj.loadDataset(nameFile);
                case 'NN'
                    filePath = fullfile("AbrilTFGfiles","Data","Lattice",p.Training,"K_NN.mat");
                    load(filePath,"K_NN");
                    filePath = fullfile("AbrilTFGfiles","Data","Lattice",p.Training,"T_NN.mat");
                    load(filePath,"T_NN","pol_deg");
            end

            RVE = cell(obj.nSubdomains(1,2),obj.nSubdomains(1,1));

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    RVE{i,j}.ndimf = 2;

                    switch p.Option
                        case 'Dataset'
                            RVE{i,j}.Kcoarse= obj.data.K{i,j};
                            RVE{i,j}.U= obj.data.T{i,j}; 
                        case 'NN'
                            RVE{i,j}.Kcoarse = computeKcoarse_NN(K_NN,[obj.tFrame(i,j),obj.tCross(i,j)],8);
                            RVE{i,j}.U       = computeT_NN_2P(obj.cellMesh{i,j},[obj.tFrame(i,j),obj.tCross(i,j)],T_NN,pol_deg);
                    end
                end
            end

            s.RVE           = RVE;
            s.mesh          = obj.createCoarseMesh();
            s.DirCond       = dir;
            s.nSubdomains   = obj.nSubdomains;
            coarseSolver    = Coarse(s);

            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.Coarsesolver = coarseSolver;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'Coarse';
            eP = Preconditioner.create(ss);
            Mcoarse = @(r) eP.apply(r);

            obj.ddDofManager=ss.ddDofManager;
            close all % Ho he afegit pq s'obren fig sense info
        end

        function NameFile=computeNameFile(obj)
            n=obj.params.nelem;
            t1=obj.tFrame;
            t2=obj.tCross;
            meshName=n+"x"+n;
            name=strings(size(t1,1),size(t1,2));

            for i=1:size(t1,1)
                for j=1:size(t1,2)
                    name(i,j) = strrep("t1_"+num2str(t1(i), '%.2f'), ".", "_")+strrep("_t2_"+num2str(t2(j), '%.2f'), ".", "_")+"-"+meshName+".mat";
                end
            end
            NameFile=name;
        end


        function loadDataset(obj,name)
            p=obj.params;
            n=p.nelem;
            Taux=cell(size(name,1),size(name,2));
            Kaux=cell(size(name,1),size(name,2));
            meshName=n+"x"+n;
            for i=1:size(name,1)
                for j=1:size(name,2)
                    filePath = fullfile("AbrilTFGfiles","Data","Lattice",p.Training,meshName,name(i,j));
                    load(filePath,"T","Kcoarse");
                    Taux{i,j}=T;
                    Kaux{i,j}=Kcoarse;
                end
            end
            obj.data.T=Taux;
            obj.data.K=Kaux;
        end


        function d = createDomainDecompositionDofManager(obj,iC,lG,bS,mR,iCR)
            s.nSubdomains     = obj.nSubdomains;
            s.interfaceConnec = iC;
            s.interfaceConnecReshaped = iCR;
            s.locGlobConnec   = lG;
            s.nBoundaryNodes  = bS{1}.mesh.nnodes;
            s.nReferenceNodes = mR.nnodes;
            s.nNodes          = obj.meshDomain.nnodes;
            s.nDimf           = obj.meshDomain.ndim;
            d = DomainDecompositionDofManager(s);
        end

        function Milu = createILUpreconditioner(~,LHS)
            s.LHS = sparse(LHS);
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);
        end

    end

    methods (Static, Access = public)

        function J = computeTotalEnergy(x,A,b)
            J = 0.5*x'*A(x)-b'*x;
        end

    end

end
