classdef TutorialEIFEM_parametric < handle

    properties (Access = public)

    end
    properties (Access = private)
        referenceMesh
        meshDomain
        boundaryConditions
        bcApplier
        LHS
        RHS
        fileNameEIFEM
        tolSameNode
        solverType
        r
    end


    properties (Access = private)
        nSubdomains
        xmin
        xmax
        ymin
        ymax
        cx
        cy
        Nr
        Ntheta
    end

    methods (Access = public)

        function obj = TutorialEIFEM_parametric()
            close all
            obj.init()

            %             obj.createReferenceMesh();
            %             bS  = obj.referenceMesh.createBoundaryMesh();
            %             [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain();
            mSbd = obj.createSubDomainMeshes();
            bS = mSbd{1,1}.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomainJoiner(mSbd);
            obj.meshDomain = mD;
            %             mD.plot()
            [bC,dir] = obj.createBoundaryConditions();
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            [LHSr,RHSr] = obj.createElasticProblem();

            LHSfun = @(x) LHSr*x;
            Meifem       = obj.createEIFEMPreconditioner(dir,iC,lG,bS,iCR,discMesh,mSbd);
            Milu         = obj.createILUpreconditioner(LHSr);
            Mmult        = @(r) Preconditioner.multiplePrec(r,LHSfun,Milu,Meifem,Milu);
            Mid          = @(r) r;


            tol = 1e-8;
            x0 = zeros(size(RHSr));
            xSol=x0;
%             xSol = LHSr\RHSr;

            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSfun,RHSr,x0,Mmult,tol,xSol,obj.meshDomain,obj.bcApplier);
            [uCG,residualCG,errCG,errAnormCG]    = PCG.solve(LHSfun,RHSr,x0,Mid,tol,xSol);
            obj.plotResidual(residualPCG,errPCG,errAnormPCG,residualCG,errCG,errAnormCG)
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [15 5]; %nx ny

            %             filePath = ['./EPFL/data_' num2str(obj.r(i), '%.3f') '.mat'];
            obj.tolSameNode = 1e-10;
            obj.solverType = 'REDUCED';
            %             obj.r  = 0.797227;
            %              obj.r  = 0.1;
            %             obj.r  = [0.1 , 0.2; 0.3, 0.4];
            rng(0);
            maxr = 0.1;
            minr = 0.1;
            obj.r= (maxr - minr) * rand(obj.nSubdomains(2),obj.nSubdomains(1)) + minr;
            obj.xmin = -1;
            obj.xmax = 1;
            obj.ymin = -1;
            obj.ymax = 1;
            obj.cx = 0;
            obj.cy = 0;
            obj.Nr=7;
            obj.Ntheta=14;
            obj.fileNameEIFEM = './EPFL/parametrizedEIFEMLagrange40.mat';
%             obj.fileNameEIFEM = './EPFL/dataEIFEMQ12_2.mat';
        end

        function createReferenceMesh(obj)
            %             obj.loadRefMesh();
            obj.referenceMesh = obj.mesh_rectangle_via_triangles(obj.r);
        end

        function  mSbd= createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            Lx = obj.xmax-obj.xmin;
            Ly = obj.ymax-obj.ymin;
            for jDom = 1:nY
                for iDom = 1:nX
                    refMesh = obj.mesh_rectangle_via_triangles(obj.r(jDom,iDom));
                    coord0 = refMesh.coord;
                    s.coord(:,1) = coord0(:,1)+Lx*(iDom-1);
                    s.coord(:,2) = coord0(:,2)+Ly*(jDom-1);
                    s.connec = refMesh.connec;
                    mIJ     = Mesh.create(s);
                    %                     plot(mIJ)
                    %                     hold on;
                    mSbd{jDom,iDom} = mIJ;
                end
            end
            obj.referenceMesh = mSbd{1,1};
        end


        function loadRefMesh(obj)
            filename = obj.fileNameEIFEM;
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            s.connec   = EIFEoper.MESH.CN;
            s.interType = 'QUADRATIC';
            obj.referenceMesh = Mesh.create(s);
        end

        function mesh = mesh_rectangle_via_triangles(obj,r)
            % Mesh rectangle with hole by meshing 4 triangles and flipping the mesh

            % Define corners
            corners = [obj.xmax obj.ymax;
                obj.xmin obj.ymax;
                obj.xmin obj.ymin;
                obj.xmax obj.ymin];

            nodes_all = [];
            elements_all = [];
            node_offset = 0;

            % Loop over 4 corners - mesh each triangle sector and transform it
            for k = 1:4
                c1 = corners(k,:);
                c2 = corners(mod(k,4)+1,:);

                % Mesh one triangle sector
                [nodes, elements] = obj.mesh_triangle_sector(obj.cx, obj.cy, c1, c2, r, obj.Nr, obj.Ntheta);

                % Append nodes, elements with offset
                elements = elements + node_offset;
                nodes_all = [nodes_all; nodes];
                elements_all = [elements_all; elements];
                node_offset = size(nodes_all,1);
            end

            tol = 1e-8;

            % Step 1: Snap coordinates to a regular grid
            rounded_nodes = round(nodes_all / tol) * tol;
            %             rounded_nodes = round(nodes_all * 1e12) / 1e12;
            [~, ia, ic] = unique(rounded_nodes, 'rows', 'stable');
            nodes_final = nodes_all(ia, :);
            elements_final = ic(elements_all);


            %             [~, ~, J] = uniquetol(nodes_all, 1e-12, 'ByRows', true);
            %             firstidx = accumarray(J, (1:length(J))', [], @min);
            %             nodes_final = nodes_all(firstidx, :);
            %             elements_final = J(elements_all);
            %
            % %             % Merge duplicate nodes (within a tolerance)
            % %             [nodes_final, ~, ic] = uniquetol(nodes_all, 1e-12, 'ByRows', true,'Stable');
            % %
            % %             % Re-map element indices to unique node list
            % %             elements_final = ic(elements_all);

            s.coord = nodes_final;
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-1e-9,-1e-9];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-1e-9,1e-9];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[1e-9,-1e-9];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[1e-9,1e-9];

            %             if isempty(obj.connec)
            %                 s.connec = elements_final;
            %                 obj.connec = s.connec;
            %             else
            %                 s.connec = obj.connec;
            %             end

            s.connec = elements_final;
            mesh = Mesh.create(s);
        end

        function [nodes, elements] = mesh_triangle_sector(obj,cx, cy, corner1, corner2, r_inner, Nr, Ntheta)
            v1 = corner1 - [cx, cy];
            v2 = corner2 - [cx, cy];

            % Prepare storage
            nodes = [];
            for i = 0:Ntheta
                % Linear interpolation between v1 and v2
                t = i / Ntheta;
                edge_vec = (1 - t) * v1 + t * v2;

                % Max radius along this direction
                r_max = norm(edge_vec);

                % Normalize direction vector
                dir = edge_vec / r_max;

                % Create points from r_inner to r_max along dir
                r_vals = linspace(r_inner, r_max, Nr + 1)';
                pts = [cx, cy] + r_vals * dir;

                nodes = [nodes; pts];
            end

            % Reshape node indices into a grid: (Nr+1) rows × (Ntheta+1) cols
            node_grid = reshape(1:size(nodes,1), Nr+1, Ntheta+1);

            % Build quad elements
            elements = [];
            for j = 1:Ntheta
                for i = 1:Nr
                    n1 = node_grid(i, j);
                    n2 = node_grid(i+1, j);
                    n3 = node_grid(i+1, j+1);
                    n4 = node_grid(i, j+1);
                    elements(end+1, :) = [n1 n2 n3 n4];
                end
            end
        end


        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj)
            s.nsubdomains   = obj.nSubdomains;
            s.meshReference = obj.referenceMesh;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE.create(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomainJoiner(obj,mSbd)
            s.nsubdomains   = obj.nSubdomains;
            s.meshReference = obj.referenceMesh;
            s.tolSameNode = obj.tolSameNode;
            s.meshSbd     = mSbd;
            %             m = MeshCreatorFromRVE.create(s);
            m = MeshJoiner(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end


        function mCoarse = createCoarseMesh(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh(mR);
%                         s.meshReference = obj.createReferenceCoarseMesh2(mR,2);
%             s.meshReference = obj.loadReferenceCoarseMesh(mR);
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE.create(s);
%             tic
            [mCoarse,~,~] = mRVECoarse.create();
%             tic
%             mCoarse2 = QuadMesh((obj.xmax-obj.xmin)*obj.nSubdomains(1),(obj.ymax-obj.ymin)*obj.nSubdomains(2),obj.nSubdomains(1),obj.nSubdomains(2));toc
        end

                function cMesh = createReferenceCoarseMesh(obj,mR)
                    xmax = max(mR.coord(:,1));
                    xmin = min(mR.coord(:,1));
                    ymax = max(mR.coord(:,2));
                    ymin = min(mR.coord(:,2));
                    coord(1,1) = xmin;  coord(1,2) = ymin;
                    coord(2,1) = xmax;  coord(2,2) = ymin;
                    coord(3,1) = xmax;  coord(3,2) = ymax;
                    coord(4,1) = xmin;  coord(4,2) = ymax;
%                     connec = [2 3 4 1];
                    connec = [1 2 3 4];
                    s.coord = coord;
                    s.connec = connec;
                    cMesh = Mesh.create(s);
                end


        function cMesh = createReferenceCoarseMesh2(obj, mR, p)
            % createReferenceCoarseMesh - Creates a reference coarse mesh with only edge nodes
            %
            % Inputs:
            %   obj - object (unused)
            %   mR  - reference mesh (contains coord)
            %   p   - order (1 = 4 nodes, 2 = 8 nodes, etc.)
            %
            % Output:
            %   cMesh - coarse mesh with edge nodes only

            % Domain bounds
            xmax = max(mR.coord(:,1));
            xmin = min(mR.coord(:,1));
            ymax = max(mR.coord(:,2));
            ymin = min(mR.coord(:,2));

            % Parametric coordinates along edges
            xi = linspace(xmin, xmax, p+1);
            eta = linspace(ymin, ymax, p+1);

            coord = [];

            % --- Bottom edge (y = ymin)
            xb = xi;
            yb = ymin * ones(1, p+1);
            coord = [coord; xb(:), yb(:)];

            % --- Right edge (x = xmax), skip first corner
            xr = xmax * ones(1, p);
            yr = eta(2:end);
            coord = [coord; xr(:), yr(:)];

            % --- Top edge (y = ymax), skip first corner (moving right→left)
            xt = xi(end-1:-1:1);
            yt = ymax * ones(1, p);
            coord = [coord; xt(:), yt(:)];

            % --- Left edge (x = xmin), skip first and last corner (moving top→bottom)
            xl = xmin * ones(1, p-1);
            yl = eta(end-1:-1:2);
            coord = [coord; xl(:), yl(:)];

            % --- Connectivity (just one polygonal element)
            connec = 1:size(coord,1);

            % Build structure
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end

         function cMesh = loadReferenceCoarseMesh(obj,mR)
            bS  = mR.createBoundaryMesh();
            bS2{1} = bS{3}; bS2{2} = bS{2}; bS2{3} = bS{4}; bS2{4} = bS{1}; % reorder boundaries
            bS = bS2;
            nbd = size(bS,2);
            interpType = [1,1,1,1];
            inode = 1;
            for ibd = 1:nbd
                maxCoord  = max(bS{ibd}.mesh.coord);
                minCoord  = min(bS{ibd}.mesh.coord);
                meanCoord = (maxCoord+minCoord)/2;
                val = ibd<=nbd/2;
                if interpType(ibd) == 1
                    coord(inode,:)   = val*minCoord + abs((val-1))*maxCoord;
                    coord(inode+1,:) = val*maxCoord + abs((val-1))*minCoord;
                    inode=inode + 2;
                else
                    coord(inode,:)   = val*minCoord + abs((val-1))*maxCoord;
                    coord(inode+1,:) = meanCoord;
                    coord(inode+2,:) = val*maxCoord + abs((val-1))*minCoord;
                    inode=inode + 3;
                end
            end

%              coord(1,:)  = [ 0.378041543026706 , -0.843442136498517 ];
%              coord(2,:)  = [ 1.49050445103858  , -0.843442136498517 ];
%              coord(3,:)  = [ 2.60296735905045  , -0.843442136498517 ];
%              coord(4,:)  = [ 2.98100890207715  ,  0                 ];
%              coord(5,:)  = [ 2.98100890207715  ,  0.314540059347181 ];
%              coord(6,:)  = [ 2.60296735905045  ,  1.1579821958457   ];
%              coord(7,:)  = [ 1.49050445103858  ,  1.1579821958457   ];
%              coord(8,:)  = [ 0.378041543026706 ,  1.1579821958457   ];
%              coord(9,:)  = [ 0                 ,  0.314540059347181 ];
%              coord(10,:) = [ 0                 ,  0                 ];
%         
%          
                     
        
            connec = [1 2 3 4 5 6 7 8];
            s.coord = coord;
            s.connec = connec;
            cMesh = Mesh.create(s);
        end



        function createBCapplier(obj)
            s.mesh                  = obj.meshDomain;
            s.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(s);
        end

        function material = createMaterial(obj,mesh)
            E  = 1;
            nu = 1/3;
            %              Epstr  = E/(1-nu^2);
            %             nupstr = nu/(1-nu);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = ConstantFunction.create(E,mesh);
            s.poisson = ConstantFunction.create(nu,mesh);
            tensor    = Material.create(s);
            material  = tensor;
        end


        function [bC,Dir] = createBoundaryConditions(obj)
            minx = min(obj.meshDomain.coord(:,1));
            maxx = max(obj.meshDomain.coord(:,1));
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;
            dirichletFun = DirichletCondition(obj.meshDomain, Dir{1});

            mesh = obj.meshDomain;
            PL.domain    = @(coor) isRight(coor);
            PL.direction = 2;
            PL.value     = -0.001;
            pointload = TractionLoad(mesh,PL,'DIRAC');
            % % %             need this because force applied in the face not in a point
            % %             pointload.values        = pointload.values/size(pointload.dofs,1);
            % %             fvalues                 = zeros(mesh.nnodes*mesh.ndim,1);
            % %             fvalues(pointload.dofs) = pointload.values;
            % %             fvalues                 = reshape(fvalues,mesh.ndim,[])';
            % %             pointload.fun.setFValues(fvalues);

            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bC             = BoundaryConditions(s);
        end



        function [LHSr,RHSr] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
            RHSr       = obj.computeForces(lhs,u);
        end

        function [LHS,LHSr] = computeStiffnessMatrix(obj,mesh,dispFun,C)
            % s.type     = 'ElasticStiffnessMatrix';
            % s.mesh     = mesh;
            % s.test     = dispFun;
            % s.trial    = dispFun;
            % s.material = mat;
            % s.quadratureOrder = 2;
            % lhs = LHSIntegrator.create(s);
            % LHS = lhs.compute();

            LHS = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),dispFun,dispFun,mesh,'Domain',2);
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
            if strcmp(obj.solverType,'REDUCED')
                bc      = obj.boundaryConditions;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    R = -stiffness(:,dirich)*dirichV;
                else
                    R = zeros(sum(obj.uFun.nDofs(:)),1);
                end
                rhs = rhs+R;
            end
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(rhs);
        end

        function RHS = computeForces1(obj,stiffness,u)
            % s.type      = 'Elastic';
            % s.scale     = 'MACRO';
            % s.dim.ndofs = u.nDofs;
            % s.BC        = obj.boundaryConditions;
            % s.mesh      = obj.meshDomain;
            % RHSint      = RHSIntegrator.create(s);
            % rhs         = RHSint.compute();
            % % Perhaps move it inside RHSint?
            % R           = RHSint.computeReactions(stiffness);
            % RHS = rhs+R;



            ndofs = u.nDofs;
            bc            = obj.boundaryConditions;
            neumann       = bc.pointload_dofs;
            neumannValues = bc.pointload_vals;
            rhs = zeros(ndofs,1);
            if ~isempty(neumann)
                rhs(neumann) = neumannValues;
            end
            if strcmp(obj.solverType,'REDUCED')
                bc      = obj.boundaryConditions;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    R = -stiffness(:,dirich)*dirichV;
                else
                    R = zeros(sum(ndofs(:)),1);
                end
                rhs = rhs+R;
            end
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(rhs);
        end

        function Meifem = createEIFEMPreconditioner(obj,dir,iC,lG,bS,iCR,dMesh,mSbd)
            mR = obj.referenceMesh;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            EIFEMfilename = obj.fileNameEIFEM;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            filename        = EIFEMfilename;
            s.RVE           = TrainedRVE(filename);

%             training 1 subdomain
                        data = Training(mR,1);
                        p = OfflineDataProcessor(data);
                        EIFEoper = p.computeROMbasis();
                         EIFEoper.Kcoarse = @(r) EIFEoper.Kcoarse;
                         EIFEoper.Udef   = @(r) EIFEoper.Udef;
                         EIFEoper.Urb   = @(r) EIFEoper.Urb;
                        s.RVE           = TrainedRVE(EIFEoper);

            % training every subdomain
            %             EIFEoper = obj.trainSubdomain(mSbd);
            %             s.RVE           = TrainedRVE(EIFEoper);

            s.mesh          = obj.createCoarseMesh(mR);
            %            s.mesh          = obj.loadCoarseMesh(mR);
            s.DirCond       = dir;
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

        function EIFEoper = trainSubdomain(obj,mSbd)
            k = 1;
            for j = 1:obj.nSubdomains(2)
                for i= 1:obj.nSubdomains(1)
                    m = mSbd{j,i};
                    data = Training(m);
                    p = OfflineDataProcessor(data);
                    EIFE = p.computeROMbasis();
                    EIFEoper.Kcoarse(:,:,k) = EIFE.Kcoarse;
                    EIFEoper.Udef(:,:,k) = EIFE.Udef;
                    EIFEoper.Urb(:,:,k) = EIFE.Urb;
                    k=k+1;
                end
            end
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
            s.LHS = LHS;
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);
        end

        function plotResidual(obj,residualPCG,errPCG,errAnormPCG,residualCG,errCG,errAnormCG)
            figure
            plot(residualPCG,'linewidth',2)
            hold on
            plot(residualCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend({'CG + ILU-EIFEM-ILU','CG'},'FontSize',12)
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

        end

    end

end
