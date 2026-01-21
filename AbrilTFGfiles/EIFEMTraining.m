classdef EIFEMTraining < handle

    properties (GetAccess = public, SetAccess = private)
        mesh
        E
        nu
        Coarseorder
    end
    properties (Access = private)
        radius
        meshDomain
        DDdofManager
        domainIndices
        Inclusion
        tolSameNode
        cellMeshes
        nSubdomains
    end


    methods (Access = public)

        function obj = EIFEMTraining(cParams)
            obj.init(cParams);
        end


        function data=train(obj)
            % create MeshDomain
            if sum(obj.nSubdomains > 1)>= 1
                obj.repeatMesh();
            else
                obj.cellMeshes= {obj.mesh};
                obj.meshDomain = obj.mesh;
            end
            
            [bMesh, lGCBd] = obj.meshDomain.createSingleBoundaryMesh();
            cF = CoarseFunction(bMesh,obj.Coarseorder);

            % ELASTIC HARMONIC EXTENSION
            s.mesh=obj.meshDomain;
            s.uFun=LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            s.lambdaFun=LagrangianFunction.create(bMesh,obj.meshDomain.ndim,'P1');
            s.material=obj.createMaterial();
            s.dirichletFun=cF.f;
            s.localGlobalConnecBd = lGCBd;
            e  = ElasticHarmonicExtension(s);
            [u,~,K] = e.solve();

            [data.uSbd,data.LHSsbd]    = obj.extractDomainData(u,K);
             data.mesh= obj.mesh;
             data.E= obj.E;
             data.nu = obj.nu;
             data.Coarseorder= obj.Coarseorder;
        end
    end

    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.radius=cParams.r;
            obj.Inclusion=cParams.Inclusion;
            switch cParams.Sampling
                case 'Isolated'
                    obj.nSubdomains   = [1 1]; %nx ny
                    obj.domainIndices = [1 1];
                case 'Oversampling'
                    obj.nSubdomains   = [3 3]; %nx ny
                    obj.domainIndices = [2 2];
            end
            obj.tolSameNode = 1e-10;
            obj.Coarseorder = 1;
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

        % function uD = AnalyticalDirCond(obj)
        %     xmax = max(obj.meshDomain.coord(:,1));
        %     ymax = max(obj.meshDomain.coord(:,2));
        %     xmin = min(obj.meshDomain.coord(:,1));
        %     ymin = min(obj.meshDomain.coord(:,2));
        %     a = (-xmin+xmax)/2;
        %     b = (-ymin+ymax)/2;
        %     x0 = xmin+a;
        %     y0 = ymin+b;
        % 
        %     f1x = @(x) [1/(4)*(1-(x(1,:,:)-x0)/a).*(1-(x(2,:,:)-y0)/b);...
        %         0*x(2,:,:)  ];
        %     f2x = @(x) [1/(4)*(1+(x(1,:,:)-x0)/a).*(1-(x(2,:,:)-y0)/b);...
        %         0*x(2,:,:)  ];
        %     f3x = @(x) [1/(4)*(1+(x(1,:,:)-x0)/a).*(1+(x(2,:,:)-y0)/b);...
        %         0*x(2,:,:)  ];
        %     f4x = @(x) [1/(4)*(1-(x(1,:,:)-x0)/a).*(1+(x(2,:,:)-y0)/b);...
        %         0*x(2,:,:)  ];
        % 
        %     f1y = @(x) [0*x(1,:,:);...
        %         1/(4)*(1-(x(1,:,:)-x0)/a).*(1-(x(2,:,:)-y0)/b) ];
        %     f2y = @(x) [0*x(1,:,:);...
        %         1/(4)*(1+(x(1,:,:)-x0)/a).*(1-(x(2,:,:)-y0)/b) ];
        %     f3y = @(x) [0*x(1,:,:);...
        %         1/(4)*(1+(x(1,:,:)-x0)/a).*(1+(x(2,:,:)-y0)/b) ];
        %     f4y = @(x) [0*x(1,:,:);...
        %         1/(4)*(1-(x(1,:,:)-x0)/a).*(1+(x(2,:,:)-y0)/b)];
        % 
        %     f     = {f1x f1y f2x f2y f3x f3y f4x f4y}; %
        %     nfun = size(f,2);
        %     for i=1:nfun
        %         uD{i}  = AnalyticalFunction.create(f{i},bMesh);
        %     end
        % end

    end


end
