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

            obj.createReferenceMesh();
            bS  = obj.referenceMesh.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain();
%             mSbd = obj.createSubDomainMeshes();
%             [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain(mSbd);
            obj.meshDomain = mD;
            mD.plot()
            [bC,dir] = obj.createBoundaryConditions();
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            [LHSr,RHSr] = obj.createElasticProblem();

            LHSfun = @(x) LHSr*x;
            Meifem       = obj.createEIFEMPreconditioner(dir,iC,lG,bS,iCR,discMesh);
            Milu         = obj.createILUpreconditioner(LHSr);
            Mmult        = @(r) Preconditioner.multiplePrec(r,LHSfun,Milu,Meifem,Milu);
            Mid          = @(r) r;


            tol = 1e-8;
            x0 = zeros(size(RHSr));
            xSol = LHSr\RHSr;

            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSfun,RHSr,x0,Mmult,tol,xSol);     
            [uCG,residualCG,errCG,errAnormCG]    = PCG.solve(LHSfun,RHSr,x0,Mid,tol,xSol);     
            obj.plotResidual(residualPCG,errPCG,errAnormPCG,residualCG,errCG,errAnormCG)
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [15 5]; %nx ny
            
%             filePath = ['/home/raul/Documents/GitHub/EPFL/data_' num2str(obj.r(i), '%.3f') '.mat'];
            obj.tolSameNode = 1e-10;
            obj.solverType = 'REDUCED';
            obj.r  = 0.8;
%             obj.r  = [0.1 , 0.2; 0.3, 0.4];
            obj.xmin = -1; 
            obj.xmax = 1;
            obj.ymin = -1;
            obj.ymax = 1;
            obj.cx = 0; 
            obj.cy = 0;
            obj.Nr=7;
            obj.Ntheta=14;
            obj.fileNameEIFEM = '/home/raul/Documents/GitHub/EPFL/parametrizedEIFEM.mat';
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
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)-[1e-9,0];
            s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)-[1e-9,0];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[1e-9,0];
            s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
                s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[1e-9,0];

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
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE.create(s);
            [mCoarse,~,~] = mRVECoarse.create();
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
            connec = [2 3 4 1];
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
            PL.value     = -0.1;
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

            LHS = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),dispFun,dispFun,mesh,2);
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);
        end

        function RHS = computeForces(obj,stiffness,u)
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

        function Meifem = createEIFEMPreconditioner(obj,dir,iC,lG,bS,iCR,dMesh)
            mR = obj.referenceMesh;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/RAUL_rve_10_may_2024/EXAMPLE/EIFE_LIBRARY/DEF_Q4porL_2s_1.mat';
            EIFEMfilename = obj.fileNameEIFEM;
            % obj.EIFEMfilename = '/home/raul/Documents/Thesis/EIFEM/05_HEXAG2D/EIFE_LIBRARY/DEF_Q4auxL_1.mat';
            filename        = EIFEMfilename;
            s.RVE           = TrainedRVE(filename);
%             data = Training(mR);
%             p = OfflineDataProcessor(data);
%             EIFEoper = p.computeROMbasis();
%             s.RVE           = TrainedRVE(EIFEoper);
            s.mesh          = obj.createCoarseMesh(mR);
%            s.mesh          = obj.loadCoarseMesh(mR);
            s.DirCond       = dir;
            s.nSubdomains = obj.nSubdomains;
            s.mu          = obj.r;
            eifem           = EIFEMnonPeriodic(s);

            ss.ddDofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            ss.EIFEMsolver = eifem;
            ss.bcApplier = obj.bcApplier;
            ss.dMesh     = dMesh;
            ss.type = 'EIFEM';
            eP = Preconditioner.create(ss);
            Meifem = @(r) eP.apply(r);
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
