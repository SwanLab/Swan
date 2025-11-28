classdef CoarseTesting_Abril< handle

    properties (Access = public)

    end
    properties (Access = private)
        nSubdomains
        rSubdomains
        mSubdomains
        ic
        icr
        lg
        bs

        meshDomain
        referenceMesh
        cellMeshes
        discMesh
        boundaryConditions
        bcApplier
        LHS
        RHS
        r
        centroids

        fileNameCorase
        tolSameNode
        NN
        data
        loadData

        xmin 
        xmax 
        ymin 
        ymax         

    end


    methods (Access = public)

        function obj = CoarseTesting_Abril()
            obj.init()
            mR = obj.createReferenceMesh();  %Crea la background mesh
            bS  = mR.createBoundaryMesh();    %crea el boundary de la mesh
            obj.referenceMesh = mR;

            if sum(obj.nSubdomains > 1)>= 1
                obj.repeatMesh();  % Crea el domini
            else
                obj.meshDomain = mR;
            end
            
            [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            [LHS,RHS,LHSf] = obj.createElasticProblem();
            obj.LHS = LHSf;
            %  LHS = 0.5*(LHS+LHS');

            LHSf = @(x) LHS*x;
            RHSf = RHS;
            Usol = LHS\RHS;
            Ufull = obj.bcApplier.reducedToFullVectorDirichlet(Usol);  %AIXO PQ SERVIEX???
            %obj.plotSolution(Ufull,obj.meshDomain,1,1,0,obj.bcApplier,0)

            Meifem       = obj.createEIFEMPreconditioner(dir,obj.ic,obj.lg,bS,obj.icr,obj.discMesh);            
            Milu         = obj.createILUpreconditioner(LHS);
%           Mcoarse      = obj.createCoarseNNPreconditioner(mR,dir,iC,lG,bS,iCR,discMesh);
            Mid            = @(r) r;

            MiluCG = @(r,iter) Preconditioner.InexactCG(r,LHSf,Milu,RHSf);

            tol = 1e-8;
            tic
            x0 = zeros(size(RHSf));

            Mmult = @(r) Preconditioner.multiplePrec(r,LHSf,Milu,Meifem,Milu);
            tic
            %           tau = @(r,A) 1;
            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol,obj.meshDomain,obj.bcApplier);
            %            [uCG,residualPCG,errPCG,errAnormPCG] = RichardsonSolver.solve(LHSf,RHSf,x0,Mmult,tol,tau,Usol);
            toc
            
            xFull = obj.bcApplier.reducedToFullVectorDirichlet(uPCG);
            s.mesh = obj.meshDomain;
            s.ndimf = obj.meshDomain.ndim;
            s.order = 'P1';
            s.fValues = reshape(xFull,2,[])';
            uFun = LagrangianFunction(s);
       
            %uFun.print('TestCoarseAbril1','Paraview');
            s.fValues = reshape(Ufull,2,[])';
            RealFun=LagrangianFunction(s);
          

            %Compute hole
            obj.computeSubdomainCentroid();
            CoarsePlotSolution(uFun, obj.meshDomain, obj.bcApplier,'TestCoarseAbril', obj.r, obj.centroids);
            CoarsePlotSolution(RealFun, obj.meshDomain, obj.bcApplier,'TestRealAbril', obj.r, obj.centroids);
            
            close all
            figure
            plot(residualPCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            xlabel('Iteration')
            ylabel('Residual')
            title("Residual PCG")

            figure
            plot(errPCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            xlabel('Iteration')
            ylabel('||error||_{L2}')
            title("error PCG")

            figure
            plot(errAnormPCG,'linewidth',2)
            hold on
            set(gca, 'YScale', 'log')
            xlabel('Iteration')
            ylabel('Energy norm')
            title("Err Anorm PCG")
        end

    end

    methods (Access = private)

        function init(obj)
            obj.loadData=true();  % true  --> DATASET
                                  % false --> NN

            obj.r= 0.1*ones(1,8); % Comentar, es per comparar amb el cas de l'Albert
            %obj.r=[0.1,0.1,0.1];

            obj.nSubdomains    = size(obj.r');
            obj.mSubdomains    = [];
            obj.tolSameNode    = 1e-10;
            
            if obj.loadData==true()  % Cargar del dataset
                nameFile=["UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat",...
                "UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat",...
                "UL_r0_1000-20x20.mat", "UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat"...
                "UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat"];
            
                %nameFile=["UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat","UL_r0_1000-20x20.mat"]
                obj.loadT(nameFile);
                obj.loadK(nameFile);
            else                    %Cargar de la NN
                nameNN= ["K_NN.mat","T_NN.mat"];
                obj.loadNN(nameNN);
            end
        end

            function loadNN(obj,nameNN)
                load(nameNN(1,1), 'K_NN');
                obj.NN.K=K_NN;
                load(nameNN(1,2), 'T_NN');
                obj.NN.T=T_NN;
            end
    
            function loadT(obj,name)
                Taux=cell(1,length(name));
                for i=1:size(name,1)
                    for j=1:size(name,2)
                        filePath = fullfile('AbrilTFGfiles', 'DataVariables', '20x20',name(i));
                        %filePath = fullfile('AlbertTFG files', 'mat files', 'Full',name(i));
                        load(filePath,"T");
                        Taux{i,j}=T;
                    end
                end
                obj.data.T=Taux;
            end
    
            function loadK(obj,name)
                Kaux=cell(1,length(name));
                for i=1:length(name)
                    filePath = fullfile('AbrilTFGfiles', 'DataVariables', '20x20',name(i));
                    %filePath = fullfile('AlbertTFG files', 'mat files', 'Full',name(i));
                    load(filePath,"K");
                    Kaux{1,i}=K;
                end
                obj.data.K=Kaux;
            end

        function repeatMesh(obj)
            [mD,mSb,iC,lG,iCR,discmesh] = obj.createMeshDomain(obj.referenceMesh);  
            obj.meshDomain = mD;        % mD:conj subdominis --> Tot el domini
            obj.cellMeshes = mSb; %??? % mSb: info de la malla a cada subdimini
            obj.ic = iC;  % info de les coordenades globals en tot el domini ???
            obj.icr = iCR; % info de les coordenades del corresponent subdomini ???
            obj.lg = lG; %??
            obj.discMesh=discmesh;
            obj.bs; 
        end

            function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj,mR)
                s.nsubdomains   = obj.nSubdomains; %nx ny
                s.meshReference = mR;
                s.tolSameNode = obj.tolSameNode;
                m = MeshCreatorFromRVE2D(s);
                [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
                close all;
            end


        function mS = createReferenceMesh(obj)
            mS = obj.createStructuredMesh();
            %lvSet    = obj.createLevelSetFunction(mS);
            %uMesh    = obj.computeUnfittedMesh(mS,lvSet);
            %mS = uMesh.createInnerMesh();

            %   mS = obj.createMeshFromGid();
            % mS = obj.createEIFEMreferenceMesh();
        end


        function mS = createMeshFromGid(obj)
            filename   = 'lattice_ex1';
            a.fileName = filename;
            femD       = FemDataContainer(a);
            mS         = femD.mesh;
        end

        function mS = createStructuredMesh(obj)
             %UnitMesh better
            x1      = linspace(-1,1,30);
            x2      = linspace(-1,1,30);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;

            obj.xmin = min(x1);            
            obj.xmax = max(x1);
            obj.ymin = min(x2);
            obj.ymax = max(x2);
            delta = 1e-9;
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-0*delta];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,0*delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[delta,-0*delta];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[delta,0*delta];


            mS = Mesh.create(s);


          
        end

        function computeSubdomainCentroid(obj)
            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                   x0=mean(obj.cellMeshes{i,j}.coord(:,1));
                   y0=mean(obj.cellMeshes{i,j}.coord(:,2));
                   obj.centroids = cat(1,obj.centroids, [x0,y0]);
                end
            end
        end

        function levelSet = createLevelSetFunction(obj,bgMesh)
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



        function mCoarse = createCoarseMesh(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh(mR);
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE2D(s);
            [mCoarse,~,~] = mRVECoarse.create();
        end


        function cMesh = createReferenceCoarseMesh(obj,mR)
            xmax = max(mR.coord(:,1)); % Obté les coordenades dels 4 vertex
            xmin = min(mR.coord(:,1));
            ymax = max(mR.coord(:,2));
            ymin = min(mR.coord(:,2));
            coord(1,1) = xmin;        % Crea els nodes i assigna als DOfs
            coord(1,2) = ymin;        % la coordenada corresponent
            coord(2,1) = xmax;
            coord(2,2) = ymin;
            coord(3,1) = xmax;
            coord(3,2) = ymax;
            coord(4,1) = xmin;
            coord(4,2) = ymax;

            connec = [1 2 3 4];    % crea conectivitats entre els 4 nodes
             connec = [ 2 3 4 1];
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
            
            material = cell(size(obj.r));

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    [young,poisson] = obj.computeElasticProperties(obj.cellMeshes{i,j}, obj.r(i,j) );

                    s.type        = 'ISOTROPIC';
                    s.ptype       = 'ELASTIC';
                    s.ndim        = obj.cellMeshes{i,j}.ndim;
                    s.young       = young;
                    s.poisson     = poisson;
                    tensor        = Material.create(s);
                    material{i,j} = tensor;

                end
            end
      

        end


        function material = createMaterialBasic(obj,mesh)
            [young,poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            tensor    = Material.create(s);
            material  = tensor;
        end

        function [young,poisson] = computeElasticProperties(obj,mesh, radius)
            E1  = 1;
            E2 = E1/1000;
            nu = 1/3;
            x0=mean(mesh.coord(:,1));
            y0=mean(mesh.coord(:,2));
%             young   = ConstantFunction.create(E,mesh);
%             poisson = ConstantFunction.create(nu,mesh);
            f   = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<radius)*E2 + ...
                        (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=radius)*E1 ; 

            young   = AnalyticalFunction.create(f,mesh);
            poisson = ConstantFunction.create(nu,mesh);
            
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

            %             Dir{2}.domain    = @(coor) isRight(coor) ;
            %             Dir{2}.direction = [2];
            %             Dir{2}.value     = 0;

%             PL.domain    = @(coor) isTop(coor);
%             PL.direction = [2];
%             PL.value     = [-0.1];
                        PL.domain    = @(coor) isRight(coor);
                        PL.direction = [2];
                        PL.value     = [1];       %Set displacement intensity ------------------------------------------------------------
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
            u = LagrangianFunction.create(obj.cellMeshes{1,1},obj.cellMeshes{1,1}.ndim,'P1');
            uBasic = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial();
            %materialBasic = obj.createMaterialBasic(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix(u,material);
            %[lhsBasic,~] = obj.computeStiffnessMatrixBasic(uBasic,materialBasic);
            RHSr       = obj.computeForces(lhs,uBasic);
        end


        function [LHS,LHSr] = computeStiffnessMatrix(obj,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.quadratureOrder = 2;
            s.test     = dispFun;
            s.trial    = dispFun;
            % LHScell    = cell(size(obj.rSubdomains));
            % LHScellGlobal = LHScell;
            LHSvect = [];

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    mesh     = obj.cellMeshes{i,j};
                    
                    C     = mat{i,j};
                    f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
                    lhs= IntegrateLHS(f,dispFun,dispFun,mesh,'Domain',2);

                    % LHScell{i,j} = full(lhs.compute());
                    LHSvect = cat(3, LHSvect, full(lhs) );
                    
                end
            end

            p.nSubdomains = obj.nSubdomains;
            p.interfaceConnec = obj.ic;
            p.interfaceConnecReshaped = obj.icr;
            p.locGlobConnec = obj.lg;
            p.nBoundaryNodes = obj.bs;
            p.nReferenceNodes = obj.cellMeshes{1,1}.nnodes;
            p.nNodes = obj.meshDomain.nnodes;
            p.nDimf = obj.meshDomain.ndim;
            

            dddm = DomainDecompositionDofManager(p);
            LHS = dddm.local2globalMatrix(LHSvect);
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);     
        end


        function [LHS,LHSr] = computeStiffnessMatrixBasic(obj,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.meshDomain;
            s.test     = dispFun;
            s.trial    = dispFun;
            s.material = mat;
            s.quadratureOrder = 2;
            lhs = LHSIntegrator.create(s);
            LHS = lhs.compute();
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

         function Meifem = createEIFEMPreconditioner(obj,dir,iC,lG,bS,iCR,dMesh,mSbd)
            mR = obj.referenceMesh;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
       %     EIFEMfilename = obj.fileNameEIFEM;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
        %    filename        = EIFEMfilename;
        %    s.RVE           = TrainedRVE(filename);
            
            %data = Training_Abril(obj,mR);
            data = Training(mR);
            p = OfflineDataProcessor(data);
            EIFEoper = p.computeROMbasis();
            s.RVE           = TrainedRVE(EIFEoper);

%             EIFEoper = obj.trainSubdomain(mSbd);
%             s.RVE           = TrainedRVE(EIFEoper);

            s.mesh          = obj.createCoarseMesh(mR);
%            s.mesh          = obj.loadCoarseMesh(mR);
            s.DirCond       = {dir};
            s.nSubdomains = obj.nSubdomains;
            s.mu          = obj.r;
            s.meshRef      = dMesh;
            eifem           = EIFEMnonPeriodic(s);

            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver = eifem;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'EIFEM';
            eP = Preconditioner.create(ss);
            Meifem = @(r) eP.apply(r);
        end        


        function Mcoarse = createCoarseNNPreconditioner(obj,mR,dir,iC,lG,bS,iCR,dMesh)
            RVE = cell(obj.nSubdomains(1,2),obj.nSubdomains(1,1));

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    RVE{i,j}.ndimf = 2;

                    if obj.loadData==true()     % Case where we load direct from dataset
                        RVE{i,j}.Kcoarse= obj.data.K{i,j};
                        RVE{i,j}.U= obj.data.T{i,j}; 
                    else                        % Case where we load from NN
                        RVE{i,j}.Kcoarse = obj.computeKcoarse(obj.r(i,j)); 
                        RVE{i,j}.U       = obj.computeTdownscaling(obj.r(i,j),obj.cellMeshes{i,j});
                    end
                end
            end

            s.RVE           = RVE;
            s.mesh          = obj.createCoarseMesh(obj.cellMeshes{1,1});
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

            close all % Ho he afegit pq s'obren fig sense info
        end


        function K=computeKcoarse(obj,r)
            K_aux1=obj.NN.K.computeOutputValues(r);
            K_aux2=zeros(8);

            idx=1;
            for n=1:8
                for m=n:8
                    K_aux2(n,m)=K_aux1(idx);
                    idx=idx+1;
                end
            end
            K=K_aux2+triu(K_aux2,1).';
        end

        function T=computeTdownscaling(obj,r,mesh)
            T=zeros(mesh.nnodes*mesh.ndim,8);
            for j=1:8                       % Constructs the 8 columns    
                Taux2=[];
                for i=1:size(mesh.coord,1)  % Evaluates all the coordenates and obtains corresponding column
                    dataInput=[r,mesh.coord(i,:)];  %
                    Taux1=obj.NN.T{1,j}.computeOutputValues(dataInput).';
                    Taux2=cat(1,Taux2,Taux1);
                end
                T(:,j)=Taux2;
            end
        end

        function Mdn = createDirichletNeumannPreconditioner(obj,mR,dir,iC,lG,bS,lhs,mSb,iCR)
            s.ddDofManager  = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            s.DirCond       = dir;
            s.bcApplier     = obj.bcApplier;
            s.LHS           = lhs;
            s.subdomainMesh = mSb;
            s.meshDomain = obj.meshDomain;
            s.type = 'DirichletNeumann';
            M = Preconditioner.create(s);
            Mdn = @(r) M.apply(r);
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

        function Milu = createILUpreconditioner(obj,LHS)
            s.LHS = sparse(LHS);
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);
        end

        function MgaussSeidel = createGaussSeidelpreconditioner(obj,LHS)
            s.LHS = LHS;
            s.type = 'GaussSeidel';
            M = Preconditioner.create(s);
            MgaussSeidel = @(r) M.apply(r);
        end

        function Mjacobi = createJacobipreconditioner(obj,LHS)
            s.LHS = LHS;
            s.type = 'Jacobi';
            M = Preconditioner.create(s);
            Mjacobi = @(r) M.apply(r);
        end

        function Mmodal = createModalpreconditioner(obj,LHS)
            s.LHS = LHS;
            s.nBasis = 8;
            s.type   = 'MODAL';
            M = Preconditioner.create(s);
            Mmodal = @(r) M.apply(r);
        end
    end

    methods (Static, Access = public)


        function J = computeTotalEnergy(x,A,b)
            J = 0.5*x'*A(x)-b'*x;
        end


        function plotSolution(x,mesh,row,col,iter,bcApplier,flag)
            if ~isempty(bcApplier)
                x = bcApplier.reducedToFullVectorDirichlet(x);
            end
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
            s.order = 'P1';
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
