classdef Training < handle

    properties (GetAccess = public, SetAccess = private)
        uSbd
        LHSsbd
        mesh
        E
        nu
        
    end
    properties (Access = private)
        meshDomain
        boundaryMeshJoined
        localGlobalConnecBd
        DirFun
        boundaryConditions
        bcApplier
        LHS
        RHS
        DDdofManager
        domainIndices

        fileNameEIFEM
        tolSameNode
        radius
    end


    properties (Access = private)
        nSubdomains
    end

    methods (Access = public)

        function obj = Training(meshRef)
            obj.init(meshRef)
            if sum(obj.nSubdomains > 1)>= 1
                obj.repeatMesh();
            else
                obj.meshDomain = obj.mesh;
            end
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.meshDomain.createSingleBoundaryMesh();
            obj.DirFun = obj.AnalyticalDirCondGeneralized(1);

            [LHS,RHS,uFun,lambdaFun] = obj.createElasticProblem();
            sol  = LHS\RHS;
            uAll = sol(1:uFun.nDofs,:);
            %             EIFEMtesting.plotSolution(full(uAll(:,1)),obj.meshDomain,1,1,1,[])
            K = LHS(1:uFun.nDofs,1:uFun.nDofs);
            [obj.uSbd,obj.LHSsbd]    = obj.extractDomainData(uAll,K);

            %             save('./Preconditioning/ROM/Training/PorousCell/OfflineData.mat','uSbd','meshRef','LHSsbd')

        end

    end

    methods (Access = private)

        function init(obj,mesh)
            obj.nSubdomains  = [5 5]; %nx ny
            obj.tolSameNode = 1e-10;
            obj.domainIndices = [3 3];
            obj.mesh = mesh;
            obj.E    = 1;
            obj.nu   = 1/3;
        end

        function repeatMesh(obj)
            bS  = obj.mesh.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain(obj.mesh);
            obj.meshDomain = mD;
            obj.DDdofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,obj.mesh,iCR);
        end

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE2D(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
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

        function [young,poisson] = computeElasticProperties(obj,mesh)
            young   = ConstantFunction.create(obj.E,mesh);
            poisson = ConstantFunction.create(obj.nu,mesh);
            %             E  = 70000;
            %             nu = 0.3;

            % for plane strain
            %             Epstr  = obj.E/(1-nu^2);
            %             nupstr = obj.nu/(1-nu);
            %             young   = ConstantFunction.create(Epstr,mesh);
            %             poisson = ConstantFunction.create(nupstr,mesh);

            % for 2 materials
            % %             E1  = 1;
            % %             E2 = E1/1000;
            % %             nu = 1/3;
            % %             x0=0;
            % %             y0=0;
            % %             young   = ConstantFunction.create(E1,mesh);
            % %             poisson = ConstantFunction.create(nu,mesh);
            % % %             f   = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<obj.radius)*E2 + ...
            % % %                         (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=obj.radius)*E1 ;
            % % % %                                      x(2,:,:).*0 ];
            % % %
            % % %             young   = AnalyticalFunction.create(f,mesh);
            % % %             poisson = ConstantFunction.create(nu,mesh);
        end

        function [LHS,RHS,u,dLambda] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            dLambda = LagrangianFunction.create(obj.boundaryMeshJoined,obj.meshDomain.ndim,'P1');
            LHS = obj.computeLHS(u,dLambda);
            RHS = obj.computeRHS(u,dLambda);
        end

        function LHS  = computeLHS(obj,u,dLambda)
            material = obj.createMaterial(obj.meshDomain);
            K = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
            C = obj.computeConditionMatrix(obj.meshDomain,u,dLambda);
            Z = zeros(dLambda.nDofs);
            LHS = [K C'; C Z];
        end

        function LHS = computeStiffnessMatrix(obj,mesh,dispFun,C)
            LHS = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),dispFun,dispFun,mesh,'Domain',2);
        end

        function C = computeConditionMatrix(obj,mesh,dispFun,dLambda)
            %             test     = LagrangianFunction.create(obj.boundaryMeshJoined, obj.meshDomain.ndim, 'P1'); % !!
            lhs = IntegrateLHS(@(u,v) DP(v,u),dLambda,dispFun,obj.meshDomain,'Boundary',2);
            C = lhs;
            %             nDofs = obj.meshDomain.nnodes*dLambda.ndimf;
            %             lhsg = sparse(nDofs,dLambda.nDofs);
            %             [iLoc,jLoc,vals] = find(lhs);
            %             l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0))';
            %             l2g_dof = l2g_dof(:);
            %             jGlob = l2g_dof(jLoc);
            %             iGlob = l2g_dof(iLoc);
            %             C = lhsg + sparse(iGlob,jLoc,vals, nDofs,dLambda.nDofs);
        end

        function RHS = computeRHS(obj,u,dLambda)
            F = zeros(u.nDofs,length(obj.DirFun));
            uD = obj.computeRHScondition(dLambda,obj.DirFun);
            RHS = [F ; uD];
        end

        function uD = computeRHScondition(obj,test,dir)
            nfun = size(dir,2);
            uD = [];
            for i=1:nfun
                d = dir{i};
                rDiri = IntegrateRHS(@(v) DP(v,d),test,obj.meshDomain,'Boundary',2);
                uD = [uD rDiri];
            end
        end

        function [u,lhs] = extractDomainData(obj,uC,LHS)
            if sum(obj.nSubdomains > 1)>= 1
                u   = obj.extractDomainDisplacements(uC);
                lhs = obj.extractDomainLHS(LHS);
            else
                u = full(uC);
                lhs = LHS;
            end
        end

        function u = extractDomainDisplacements(obj,uC)
            ntest = size(uC,2);
            for i = 1:ntest
                uD    = obj.DDdofManager.global2local(uC(:,i));
                ind   = (obj.domainIndices(1)-1)*obj.nSubdomains(1)+obj.domainIndices(2);
                u(:,i)= uD(:,ind);
            end
        end

        function lhs = extractDomainLHS(obj,LHS)
            lhs = obj.DDdofManager.global2localMatrix(LHS);
            ind = (obj.domainIndices(1)-1)*obj.nSubdomains(1)+obj.domainIndices(2);
            lhs = lhs(:,:,ind);
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

        function uD = AnalyticalDirCond(obj)
            xmax = max(obj.meshDomain.coord(:,1));
            ymax = max(obj.meshDomain.coord(:,2));
            xmin = min(obj.meshDomain.coord(:,1));
            ymin = min(obj.meshDomain.coord(:,2));
            a = (-xmin+xmax)/2;
            b = (-ymin+ymax)/2;
            x0 = xmin+a;
            y0 = ymin+b;

            %             f1x = @(x) [1/(4)*(1-x(1,:,:)).*(1-x(2,:,:));...
            %                     0*x(2,:,:)  ];
            %             f2x = @(x) [1/(4)*(1+x(1,:,:)).*(1-x(2,:,:));...
            %                     0*x(2,:,:)  ];
            %             f3x = @(x) [1/(4)*(1+x(1,:,:)).*(1+x(2,:,:));...
            %                     0*x(2,:,:)  ];
            %             f4x = @(x) [1/(4)*(1-x(1,:,:)).*(1+x(2,:,:));...
            %                     0*x(2,:,:)  ];
            %
            %             f1y = @(x) [0*x(1,:,:);...
            %                     1/(4)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
            %             f2y = @(x) [0*x(1,:,:);...
            %                     1/(4)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
            %             f3y = @(x) [0*x(1,:,:);...
            %                     1/(4)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
            %             f4y = @(x) [0*x(1,:,:);...
            %                     1/(4)*(1-x(1,:,:)).*(1+x(2,:,:))  ];

            f1x = @(x) [1/(4)*(1-(x(1,:,:)-x0)/a).*(1-(x(2,:,:)-y0)/b);...
                0*x(2,:,:)  ];
            f2x = @(x) [1/(4)*(1+(x(1,:,:)-x0)/a).*(1-(x(2,:,:)-y0)/b);...
                0*x(2,:,:)  ];
            f3x = @(x) [1/(4)*(1+(x(1,:,:)-x0)/a).*(1+(x(2,:,:)-y0)/b);...
                0*x(2,:,:)  ];
            f4x = @(x) [1/(4)*(1-(x(1,:,:)-x0)/a).*(1+(x(2,:,:)-y0)/b);...
                0*x(2,:,:)  ];

            f1y = @(x) [0*x(1,:,:);...
                1/(4)*(1-(x(1,:,:)-x0)/a).*(1-(x(2,:,:)-y0)/b) ];
            f2y = @(x) [0*x(1,:,:);...
                1/(4)*(1+(x(1,:,:)-x0)/a).*(1-(x(2,:,:)-y0)/b) ];
            f3y = @(x) [0*x(1,:,:);...
                1/(4)*(1+(x(1,:,:)-x0)/a).*(1+(x(2,:,:)-y0)/b) ];
            f4y = @(x) [0*x(1,:,:);...
                1/(4)*(1-(x(1,:,:)-x0)/a).*(1+(x(2,:,:)-y0)/b)];

            f     = {f1x f1y f2x f2y f3x f3y f4x f4y}; %
            nfun = size(f,2);
            for i=1:nfun
                uD{i}  = AnalyticalFunction.create(f{i},obj.boundaryMeshJoined);
            end
        end

        function uD = AnalyticalDirCondGeneralized(obj, order)
            % order: polynomial order (1 = linear, 2 = quadratic, 3 = cubic)
            if nargin < 2
                order = 1;
            end

            % Mesh normalization
            xmax = max(obj.meshDomain.coord(:,1));
            ymax = max(obj.meshDomain.coord(:,2));
            xmin = min(obj.meshDomain.coord(:,1));
            ymin = min(obj.meshDomain.coord(:,2));
            a = (xmax - xmin)/2;
            b = (ymax - ymin)/2;
            x0 = xmin + a;
            y0 = ymin + b;

            % 1D reference node positions
            xi = linspace(-1, 1, order+1);
            L = cell(order+1, 1);
            for i = 1:order+1
                L{i} = @(s) obj.lagrangeBasis(s, xi, i);
            end

            % Loop only over boundary node indices
            f = {};
            for i = 1:order+1
                for j = 1:order+1
                    if obj.isBoundaryNode(i, j, order+1)
                        N = @(x) L{i}((x(1,:,:)-x0)/a) .* L{j}((x(2,:,:)-y0)/b);

                        % x-component
                        fx = @(x) [N(x); 0*x(2,:,:)];

                        % y-component
                        fy = @(x) [0*x(1,:,:); N(x)];

                        f{end+1} = fx;
                        f{end+1} = fy;
                    end
                end
            end

            % Create AnalyticalFunction objects
            nfun = length(f);
            for k = 1:nfun
                uD{k} = AnalyticalFunction.create(f{k}, obj.boundaryMeshJoined);
                %                 uD{k}.plot()
            end
        end

        % Helper: Lagrange basis
        function val = lagrangeBasis(obj,s, nodes, i)
            n = length(nodes);
            val = ones(size(s));
            for j = 1:n
                if j ~= i
                    val = val .* (s - nodes(j)) / (nodes(i) - nodes(j));
                end
            end
        end

        % Helper: check if node index (i,j) is on boundary
        function isB = isBoundaryNode(obj,i, j, n)
            % n = order + 1
            isB = (i == 1 || i == n || j == 1 || j == n);
        end



    end

end
