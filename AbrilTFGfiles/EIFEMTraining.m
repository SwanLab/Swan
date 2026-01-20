classdef EIFEMTraining < handle

    properties (GetAccess = public, SetAccess = private)
        uSbd
        LHSsbd
        mesh
        E
        nu
        Coarseorder
    end
    properties (Access = private)
        radius
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
        Inclusion
        material

        fileNameEIFEM
        tolSameNode
        cellMeshes
        nSubdomains
    end


    methods (Access = public)

        function obj = EIFEMTraining(meshRef,r,params)
            obj.init(meshRef,r,params)

            if sum(obj.nSubdomains > 1)>= 1
                obj.repeatMesh();
            else
                obj.cellMeshes= {obj.mesh};
                obj.meshDomain = obj.mesh;
            end
            
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.meshDomain.createSingleBoundaryMesh();
            cF = CoarseFunction(obj.boundaryMeshJoined,obj.Coarseorder);
            obj.DirFun = cF.f;
            
            obj.material=obj.createMaterial();
            s.mesh=obj.meshDomain;
            [LHS,RHS,uFun,~] = obj.createElasticProblem();
            sol  = LHS\RHS;
            uAll = sol(1:uFun.nDofs,:);
            K = LHS(1:uFun.nDofs,1:uFun.nDofs);
            [obj.uSbd,obj.LHSsbd]    = obj.extractDomainData(uAll,K);

            % outputs
            % data.mesh= obj.mesh;
            % data.uSbd= obj.uSbd;
            % data.LHSsbd= obj.LHSsbd;
            % data.E= obj.E;
            % data.nu = obj.nu;
            % data.Coarseorder= obj.Coarseorder;

        end

    end

    methods (Access = private)

        function init(obj,mesh,r,params)
            switch params.Sampling
                case 'Isolated'
                    obj.nSubdomains   = [1 1]; %nx ny
                    obj.domainIndices = [1 1];
                case 'Oversampling'
                    obj.nSubdomains   = [3 3]; %nx ny
                    obj.domainIndices = [2 2];
            end
            
            obj.radius=r;
            obj.tolSameNode = 1e-10;
            obj.mesh = mesh;
            obj.Coarseorder = 1;
            obj.Inclusion=params.Inclusion;
            
        end

        function repeatMesh(obj)
            bS  = obj.mesh.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,~] = obj.createMeshDomain(obj.mesh);
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

        function material = createMaterial(obj)
                    [young,poisson] = obj.computeElasticProperties();
                    s.type          = 'ISOTROPIC';
                    s.ptype         = 'ELASTIC';
                    s.ndim          = obj.mesh.ndim;
                    s.young         = young;
                    s.poisson       = poisson;
                    tensor          = Material.create(s);
                    material        = tensor;
        end


        function [young,poisson] = computeElasticProperties(obj)
            obj.E  = 1;
            obj.nu  = 1/3;
            r   = obj.radius;

            switch obj.Inclusion
                case {'Hole','HoleRaul'}
                    young   = ConstantFunction.create(obj.E,obj.meshDomain);
                    poisson = ConstantFunction.create(obj.nu,onj.meshDomain);
                case 'Material'
                    E2 = obj.E/1000;
                    xmax = max(obj.mesh.coord(:,1));
                    ymax = max(obj.mesh.coord(:,2));
                    xmin = min(obj.mesh.coord(:,1));
                    ymin = min(obj.mesh.coord(:,2));
                    Lx = xmax-xmin;
                    Ly = ymax-ymin;

                    f = @(x) ...
                        ( sqrt( ...
                        (mod(x(1,:,:) - xmin, Lx) - Lx/2).^2 + ...
                        (mod(x(2,:,:) - ymin, Ly) - Ly/2).^2 ) < r ) * E2 + ...
                        ( sqrt( ...
                        (mod(x(1,:,:) - xmin, Lx) - Lx/2).^2 + ...
                        (mod(x(2,:,:) - ymin, Ly) - Ly/2).^2 ) >= r ) * obj.E;
                
                    young   = AnalyticalFunction.create(f, obj.meshDomain);
                    poisson = ConstantFunction.create(obj.nu, obj.meshDomain);
            end   
        end


        function [LHS,RHS,u,dLambda] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            dLambda = LagrangianFunction.create(obj.boundaryMeshJoined,obj.meshDomain.ndim,'P1');
            LHS = obj.computeLHS(u,dLambda);
            RHS = obj.computeRHS(u,dLambda);
        end

        function LHS  = computeLHS(obj,u,dLambda)
            material = obj.createMaterial();
            K = obj.computeStiffnessMatrix(u,material);
            C = obj.computeConditionMatrix(u,dLambda);
            Z = zeros(dLambda.nDofs);
            LHS = [K C'; C Z];
        end


        function LHS= computeStiffnessMatrix(obj,dispFun,C)
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            LHS= IntegrateLHS(f,dispFun,dispFun,obj.meshDomain,'Domain',2);
        end


        function C = computeConditionMatrix(obj,dispFun,dLambda)
            lhs = IntegrateLHS(@(u,v) DP(v,u),dLambda,dispFun,obj.meshDomain,'Boundary',2);
            C = lhs;
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

        function [centroids] =computeCentroid(obj)
            x0=mean(obj.mesh.coord(:,1));
            y0=mean(obj.mesh.coord(:,2));
            centroids = [x0,y0];
        end

    end


end
