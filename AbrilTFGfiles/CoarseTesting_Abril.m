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

        tolSameNode
        NN
        data
        params

        xmin 
        xmax 
        ymin 
        ymax         

    end


    methods (Access = public)

        function obj = CoarseTesting_Abril()
            obj.init()

            % Mesh Raul: 
            %prova=load('meshRaul.mat','refMesh')
            %mR=prova.refMesh;
            %obj.xmin =-1;
            %obj.xmax = 1;
            %obj.ymin =-1;
            %obj.ymax = 1;
    

            % COMPUTE MESH  
            mR  = obj.createReferenceMesh();  % Crea la reference mesh
            bS  = mR.createBoundaryMesh();   % Crea el boundary de la mesh
            obj.referenceMesh = mR;          % Guarda la reference Mesh   
            obj.repeatMesh();                % Crea el domini
            
            [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            [LHS,RHS,LHSf] = obj.createElasticProblem();
            obj.LHS        = LHSf;


            % EXACT SOLUTION
            LHSf   = @(x) LHS*x;
            RHSf   = RHS;
            Usol   = LHS\RHS;
            Ufull  = obj.bcApplier.reducedToFullVectorDirichlet(Usol); 


            % PRECONDITIONERS
            Meifem       = obj.createEIFEMPreconditioner(dir,obj.ic,obj.lg,bS,obj.icr,obj.discMesh);            
            Milu         = obj.createILUpreconditioner(LHS);
            Mcoarse     = obj.createCoarseNNPreconditioner(mR,dir,obj.ic,obj.lg,bS,obj.icr,obj.discMesh);
            %MiluCG      = @(r,iter) Preconditioner.InexactCG(r,LHSf,Milu,RHSf);
            Mmult        = @(r) Preconditioner.multiplePrec(r,LHSf,Milu,Meifem,Milu);


            % SOLVE THE CASE
            tol = 1e-8;
            x0  = zeros(size(RHSf));
            tic
            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x0,Mmult,tol,Usol,obj.meshDomain,obj.bcApplier);
            %   tau = @(r,A) 1;
            %   [uCG,residualPCG,errPCG,errAnormPCG] = RichardsonSolver.solve(LHSf,RHSf,x0,Mmult,tol,tau,Usol);
            toc
            xFull = obj.bcApplier.reducedToFullVectorDirichlet(uPCG);


            % EXPORT TO PARAVIEW
            s.mesh = obj.meshDomain;
            s.ndimf = obj.meshDomain.ndim;
            s.order = 'P1';
            s.fValues = reshape(xFull,2,[])';
            uFun = LagrangianFunction(s);
            
            uFun.print('ProvaIter1','Paraview');

            s.fValues = reshape(Ufull,2,[])';
            RealFun=LagrangianFunction(s);

            %obj.computeSubdomainCentroid();
            %CoarsePlotSolution(uFun, obj.meshDomain, obj.bcApplier,'TestCoarseAbril', obj.r, obj.centroids);
            %CoarsePlotSolution(RealFun, obj.meshDomain, obj.bcApplier,'TestRealAbril', obj.r, obj.centroids);
            

            % PLOTS
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        function init(obj)
            % Case Parameters
            p.Inclusion = 'Hole';         % 'Hole'/'Material'
            p.Sampling  = 'Oversampling'; % 'Isolated'/'Oversampling'
            p.loadData  = 'Dataset';      % 'Dataset'/'NN'
            p.nelem     =  20;            %  Mesh refining
            obj.params  =  p;

            % Definition of Subdomain
            obj.r              = ones(5,15)*0.3;
            obj.nSubdomains    = size(obj.r');
            obj.mSubdomains    = [];
            obj.tolSameNode    = 1e-10;
            
            % Load the data of the case
            switch p.loadData
                case 'Dataset'
                    nameFile=obj.computeNameFile();
                    obj.loadT(nameFile,p.Inclusion);
                    obj.loadK(nameFile,p.Inclusion);  
                case 'NN'
                    nameNN= ["K_NN.mat","T_NN.mat"];
                    obj.loadNN(nameNN);
            end
        end

        function NameFile=computeNameFile(obj)
            n=obj.params.nelem;
            rad=obj.r;
            meshName=n+"x"+n;
            name=strings(size(rad,1),size(rad,2));
            for i=1:size(rad,1)
                for j=1:size(rad,2)
                    name(i,j) = strrep("UL_r"+num2str(rad(i,j), '%.4f'), ".", "_")+"-"+meshName+".mat";
                end
            end
            NameFile=name;
        end

        function loadNN(obj,nameNN)
            load(nameNN(1,1), 'K_NN');
            obj.NN.K=K_NN;
            load(nameNN(1,2), 'T_NN');
            obj.NN.T=T_NN;
        end

        function loadT(obj,name,inclusion)
            n=obj.params.nelem;
            Taux=cell(size(name,1),size(name,2));
            meshName=n+"x"+n;
            for i=1:size(name,1)
                for j=1:size(name,2)
                    switch inclusion
                        case 'Material'
                                filePath = fullfile('AbrilTFGfiles', 'Data',meshName,name(i,j));
                        case 'Hole'
                                filePath = fullfile('AbrilTFGfiles', 'Data','hole',name(i,j));
                    end
                    load(filePath,"T");
                    Taux{i,j}=T;
                end
            end
            obj.data.T=Taux;
        end

        function loadK(obj,name,inclusion)
            n=obj.params.nelem;
            Kaux=cell(1,length(name));
            meshName=n+"x"+n;
            for i=1:size(name,1)
                for j=1:size(name,2)
                    switch inclusion
                        case 'Material'
                                filePath = fullfile('AbrilTFGfiles', 'Data',meshName,name(i,j));
                        case 'Hole'
                                filePath = fullfile('AbrilTFGfiles', 'Data','hole',name(i,j));
                    end
                    load(filePath,"K");
                    Kaux{i,j}=K;
                end
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
            p=obj.params;
            switch p.Inclusion
                case 'Material'
                    mS = obj.createStructuredMesh();
                case 'Hole'
                    mS = obj.createStructuredMesh();
                    lvSet    = obj.createLevelSetFunction(mS);
                    uMesh    = obj.computeUnfittedMesh(mS,lvSet);
                    mS = uMesh.createInnerMesh();
            end
        end


        function mS = createStructuredMesh(obj)
            n =obj.params.nelem;

            x1      = linspace(-1,1,n);
            x2      = linspace(-1,1,n);
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
            sLS.radius      = obj.r(1,1);
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



        function [young,poisson] = computeElasticProperties(obj,mesh, radius)
            E1  = 1;
            nu = 1/3;
            p=obj.params;

            switch p.Inclusion
                case 'Hole'
                    young   = ConstantFunction.create(E1,mesh);
                    poisson = ConstantFunction.create(nu,mesh);
                case 'Material'
                    E2 = E1/1000;
                    x0=mean(mesh.coord(:,1));
                    y0=mean(mesh.coord(:,2));
                    f   = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<radius)*E2 + ...
                                (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=radius)*E1 ; 
                    
                    young   = AnalyticalFunction.create(f,mesh);
                    poisson = ConstantFunction.create(nu,mesh);
            end            
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
            [lhs,LHSr] = obj.computeStiffnessMatrix(u,material);
            RHSr       = obj.computeForces(lhs,uBasic);
        end


        function [LHS,LHSr] = computeStiffnessMatrix(obj,dispFun,mat)
            LHSvect = [];
            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    mesh     = obj.cellMeshes{i,j};  
                    C     = mat{i,j};
                    f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
                    lhs= IntegrateLHS(f,dispFun,dispFun,mesh,'Domain',2);
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

         function Meifem = createEIFEMPreconditioner(obj,dir,iC,lG,bS,iCR,dMesh)
            mR = obj.referenceMesh;
            Data = OversamplingTraining(mR,obj.params);
            p = OfflineDataProcessor(Data);
            EIFEoper = p.computeROMbasis();
            s.RVE           = TrainedRVE(EIFEoper);
            s.mesh          = obj.createCoarseMesh();
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
            p=obj.params;
            RVE = cell(obj.nSubdomains(1,2),obj.nSubdomains(1,1));

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    RVE{i,j}.ndimf = 2;

                    switch p.loadData
                        case 'Dataset'
                            RVE{i,j}.Kcoarse= obj.data.K{i,j};
                            RVE{i,j}.U= obj.data.T{i,j}; 
                        case 'NN'
                            RVE{i,j}.Kcoarse = obj.computeKcoarse(obj.r(i,j)); 
                            RVE{i,j}.U       = obj.computeTdownscaling(obj.r(i,j),obj.cellMeshes{i,j});
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
                    dataInput=[r,mesh.coord(i,:)];  
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
