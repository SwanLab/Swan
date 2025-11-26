classdef Training_Abril < handle

    properties (GetAccess = public, SetAccess = private)
        uSbd
        LHSsbd
        mesh
        E
        nu
        Coarseorder
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
        cellMeshes
    end


    properties (Access = private)
        nSubdomains
    end

    methods (Access = public)

        function obj = Training_Abril(s,meshRef)
            obj.init(s,meshRef)
            if sum(obj.nSubdomains > 1)>= 1
                obj.repeatMesh();
            else
                obj.meshDomain = obj.mesh;
            end
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.meshDomain.createSingleBoundaryMesh();
            cF = CoarseFunction(obj.boundaryMeshJoined,obj.Coarseorder);
            obj.DirFun = cF.f;
%             obj.DirFun = obj.AnalyticalDirCondGeneralized(1);

            [LHS,RHS,uFun,lambdaFun] = obj.createElasticProblem();
            sol  = LHS\RHS;
            uAll = sol(1:uFun.nDofs,:);
                        EIFEMtesting.plotSolution(full(uAll(:,1)),obj.meshDomain,1,1,1,[])
            K = LHS(1:uFun.nDofs,1:uFun.nDofs);
            [obj.uSbd,obj.LHSsbd]    = obj.extractDomainData(uAll,K);

            %             save('./Preconditioning/ROM/Training/PorousCell/OfflineData.mat','uSbd','meshRef','LHSsbd')

        end

    end

    methods (Access = private)

        function init(obj,s,mesh)
            obj.nSubdomains  = s.nSubdomains; %nx ny
            obj.tolSameNode = 1e-10;
            obj.domainIndices = [3 3];
            obj.mesh = mesh;
            obj.E    = 1;
            obj.nu   = 1/3;
            obj.Coarseorder = 1;
        end

        function repeatMesh(obj)
            bS  = obj.mesh.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain(obj.mesh);
            obj.cellMeshes = mSb;
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

%         function material = createMaterial(obj,mesh)
%             [young,poisson] = obj.computeElasticProperties(mesh);
%             s.type    = 'ISOTROPIC';
%             s.ptype   = 'ELASTIC';
%             s.ndim    = mesh.ndim;
%             s.young   = young;
%             s.poisson = poisson;
%             tensor    = Material.create(s);
%             material  = tensor;
%         end

        function material = createMaterial(obj)

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    [young,poisson] = obj.computeElasticProperties(obj.cellMeshes{i,j} );

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


        function [young,poisson] = computeElasticProperties(obj,mesh)
%             young   = ConstantFunction.create(obj.E,mesh);
%             poisson = ConstantFunction.create(obj.nu,mesh);
            %             E  = 70000;
            %             nu = 0.3;

            % for plane strain
            %             Epstr  = obj.E/(1-nu^2);
            %             nupstr = obj.nu/(1-nu);
            %             young   = ConstantFunction.create(Epstr,mesh);
            %             poisson = ConstantFunction.create(nupstr,mesh);

            E1  = 1;
            E2 = E1/1000;
            nu = 1/3;
           radius = 0.1;
           x0=mean(mesh.coord(:,1));
            y0=mean(mesh.coord(:,2));
%             young   = ConstantFunction.create(E,mesh);
%             poisson = ConstantFunction.create(nu,mesh);
            f   = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<radius)*E2 + ...
                        (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=radius)*E1 ; 

            young   = AnalyticalFunction.create(f,mesh);
            poisson = ConstantFunction.create(nu,mesh);
        end

        function [LHS,RHS,uGlobal,dLambda] = createElasticProblem(obj)
            uGlobal = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            uLocal = LagrangianFunction.create(obj.cellMeshes{1,1},obj.cellMeshes{1,1}.ndim,'P1');
            dLambda = LagrangianFunction.create(obj.boundaryMeshJoined,obj.meshDomain.ndim,'P1');
            LHS = obj.computeLHS(uLocal,uGlobal,dLambda);
            RHS = obj.computeRHS(uGlobal,dLambda);
        end

        function LHS  = computeLHS(obj,uLocal,uGlobal,dLambda)
            material = obj.createMaterial();
            K = obj.computeStiffnessMatrix(uLocal,material);
            C = obj.computeConditionMatrix(obj.meshDomain,uGlobal,dLambda);
            Z = zeros(dLambda.nDofs);
            LHS = [K C'; C Z];
        end

%         function LHS = computeStiffnessMatrix(obj,mesh,dispFun,C)
%             LHS = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),dispFun,dispFun,mesh,'Domain',2);
%         end

        function [LHS,LHSr] = computeStiffnessMatrix(obj,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.quadratureOrder = 2;
            s.test     = dispFun;
            s.trial    = dispFun;
            % LHScell    = cell(size(obj.rSubdomains));
            % LHScellGlobal = LHScell;
            LHSl = [];

            for i = 1:obj.nSubdomains(1,2)
                for j = 1:obj.nSubdomains(1,1)
                    mesh     = obj.cellMeshes{i,j};
                    
                    C     = mat{i,j};
                    f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
                    lhs= IntegrateLHS(f,dispFun,dispFun,mesh,'Domain',2);

                    % LHScell{i,j} = full(lhs.compute());
                    LHSl = cat(3, LHSl, full(lhs) );
                    
                end
            end
            LHS = obj.DDdofManager.local2globalMatrix(LHSl);
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

    end

end
