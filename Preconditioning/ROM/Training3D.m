classdef Training3D < handle

    properties (Access = public)

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

    end


    properties (Access = private)
        nSubdomains
    end

    methods (Access = public)

        function obj = Training3D()
            close all
            obj.init()

            mR = obj.createReferenceMesh();
            bS  = mR.createBoundaryMesh();
            [mD,mSb,iC,lG,iCR,discMesh] = obj.createMeshDomain(mR);
            obj.meshDomain = mD;
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.meshDomain.createSingleBoundaryMesh3D();
            obj.DDdofManager = obj.createDomainDecompositionDofManager(iC,lG,bS,mR,iCR);
            obj.DirFun = obj.AnalyticalDirCond();


            [LHS,RHS,uFun,lambdaFun] = obj.createElasticProblem();
            sol  = LHS\RHS;
            uAll = sol(1:uFun.nDofs,:);
            [uSbd,LHSsbd]    = obj.extractDomainData(uAll,LHS);

            save('./Preconditioning/ROM/Training/PorousCell/OfflineData.mat','uSbd','mR','LHSsbd')

        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [2 2 2]; %nx ny
            %             obj.fileNameEIFEM = 'DEF_Q4auxL_1.mat';
            obj.fileNameEIFEM = 'DEF_Q4porL_1.mat';
            obj.tolSameNode = 1e-6;
            obj.domainIndices = [3 3 3];
        end

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE3D(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end

        function mS = createReferenceMesh(obj)
                mS = obj.createStructuredMesh();
%             mS = obj.createMeshFromGid();
            %             mS = obj.createEIFEMreferenceMesh();
        end

        function mS = createMeshFromGid(obj)
            %             filename   = 'lattice_ex1';
            %             a.fileName = filename;
            %             femD       = FemDataContainer(a);
            %             mS         = femD.mesh;

            load('por3Dmesh.mat');
            s.connec = porMesh.connec ;
            s.coord  = porMesh.coord;
            maxC= max(s.coord);
            minC = min(s.coord);
            shift = 1e-5;
            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,shift];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,shift];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,shift];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,shift];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,shift];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,shift];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,shift];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,shift];

            mS = Mesh.create(s);
        end

        function mS = createStructuredMesh(obj)
% % %             % Generate coordinates
% % %             x1 = linspace(0,1,2);
% % %             x2 = linspace(0,1,2);
% % %             % Create the grid
% % %             [xv,yv] = meshgrid(x1,x2);
% % %             % Triangulate the mesh to obtain coordinates and connectivities
% % %             [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');
% % % 
% % %             s.coord    = coord(:,1:2);
% % %             s.connec   = F;
% % %             mS         = Mesh.create(s);
            mS = TetraMesh(1,1,1,2,2,2);
            s.coord = mS.coord;
            s.connec = mS.connec;
             maxC= max(s.coord);
            minC = min(s.coord);
           shift = 1e-2;
            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,shift];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,shift];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,shift];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,shift];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,shift];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,shift];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,shift];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,shift];


            mS = Mesh.create(s);
        end

        function mS = createEIFEMreferenceMesh(obj)
            filename = obj.fileNameEIFEM;
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            isMin = s.coord==min(s.coord);
            isMax = s.coord==max(s.coord);

            s.connec   = EIFEoper.MESH.CN;
            mS         = Mesh.create(s);
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
            E  = 1;
            nu = 1/3;
            Epstr  = E/(1-nu^2);
            nupstr = nu/(1-nu);
            young   = ConstantFunction.create(Epstr,mesh);
            poisson = ConstantFunction.create(nupstr,mesh);
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

            C = obj.computeConditionMatrix(dLambda);
            Z = zeros(size(C,2));
            LHS = [K C; C' Z];
        end

        function LHS = computeStiffnessMatrix(obj,mesh,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.test     = dispFun;
            s.trial    = dispFun;
            s.material = mat;
            s.quadratureOrder = 2;
            lhs = LHSIntegrator.create(s);
            LHS = lhs.compute();
        end

        function C = computeConditionMatrix(obj,dLambda)
            s.type                  = 'ConditionMatrix';
            s.quadType              = 2;
            s.boundaryMeshJoined    = obj.boundaryMeshJoined;
            s.localGlobalConnecBd   = obj.localGlobalConnecBd;
            s.nnodes                = obj.meshDomain.nnodes;
            lhs      = LHSIntegrator_condition_shape_shape(s);
            test     = LagrangianFunction.create(obj.boundaryMeshJoined, obj.meshDomain.ndim, 'P1'); % !!
            C        = lhs.compute(dLambda,test);
        end

        function RHS = computeRHS(obj,u,dLambda)
            F = zeros(u.nDofs,length(obj.DirFun));
            uD = obj.computeRHScondition(dLambda,obj.DirFun);
            RHS = [F ; uD];
        end

        function uD = computeRHScondition(obj,test,dir)
            s.mesh = obj.boundaryMeshJoined;
            s.quadType = 2;
            rhs = RHSIntegratorShapeFunction(s);
            nfun = size(dir,2);
            uD = [];
            for i=1:nfun
                rDiri = rhs.compute(dir{i},test);
                uD = [uD rDiri];
            end
        end

        function [u,lhs] = extractDomainData(obj,uC,LHS)
            u   = obj.extractDomainDisplacements(uC);
            lhs = obj.extractDomainLHS(LHS);
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

        %         function uD = AnalyticalDirCond(obj)
        %             % Get bounding box
        %             xmax = max(obj.meshDomain.coord(:,1));
        %             ymax = max(obj.meshDomain.coord(:,2));
        %             zmax = max(obj.meshDomain.coord(:,3));
        %             xmin = min(obj.meshDomain.coord(:,1));
        %             ymin = min(obj.meshDomain.coord(:,2));
        %             zmin = min(obj.meshDomain.coord(:,3));
        %
        %             % Center and scaling
        %             a = (xmax - xmin)/2;
        %             b = (ymax - ymin)/2;
        %             c = (zmax - zmin)/2;
        %             x0 = xmin + a;
        %             y0 = ymin + b;
        %             z0 = zmin + c;
        %
        %             % Shape function gradients (as vector fields)
        %             % Each function returns a 3xN matrix (vector field)
        %             f1x = @(x) [ (1 - (x(1,:,:)-x0)/a - (x(2,:,:)-y0)/b - (x(3,:,:)-z0)/c); 0*x(1,:,:); 0*x(1,:,:) ];
        %             f1y = @(x) [ 0*x(1,:,:); (1 - (x(1,:,:)-x0)/a - (x(2,:,:)-y0)/b - (x(3,:,:)-z0)/c); 0*x(1,:,:) ];
        %             f1z = @(x) [ 0*x(1,:,:); 0*x(1,:,:); (1 - (x(1,:,:)-x0)/a - (x(2,:,:)-y0)/b - (x(3,:,:)-z0)/c) ];
        %
        %             f2x = @(x) [ (x(1,:,:)-x0)/a; 0*x(1,:,:); 0*x(1,:,:) ];
        %             f2y = @(x) [ 0*x(1,:,:); (x(1,:,:)-x0)/a; 0*x(1,:,:) ];
        %             f2z = @(x) [ 0*x(1,:,:); 0*x(1,:,:); (x(1,:,:)-x0)/a ];
        %
        %             f3x = @(x) [ (x(2,:,:)-y0)/b; 0*x(1,:,:); 0*x(1,:,:) ];
        %             f3y = @(x) [ 0*x(1,:,:); (x(2,:,:)-y0)/b; 0*x(1,:,:) ];
        %             f3z = @(x) [ 0*x(1,:,:); 0*x(1,:,:); (x(2,:,:)-y0)/b ];
        %
        %             f4x = @(x) [ (x(3,:,:)-z0)/c; 0*x(1,:,:); 0*x(1,:,:) ];
        %             f4y = @(x) [ 0*x(1,:,:); (x(3,:,:)-z0)/c; 0*x(1,:,:) ];
        %             f4z = @(x) [ 0*x(1,:,:); 0*x(1,:,:); (x(3,:,:)-z0)/c ];
        %
        %             % Wrap into array of functions
        %             f = {f1x, f1y, f1z, ...
        %                 f2x, f2y, f2z, ...
        %                 f3x, f3y, f3z, ...
        %                 f4x, f4y, f4z};
        %
        %             % Generate AnalyticalFunction objects
        %             nfun = numel(f);
        %             for i = 1:nfun
        %                 uD{i} = AnalyticalFunction.create(f{i}, obj.meshDomain.ndim, obj.boundaryMeshJoined);
        %             end
        %         end

        function uD = AnalyticalDirCond(obj)
            % Get bounding box
            xmax = max(obj.meshDomain.coord(:,1));
            ymax = max(obj.meshDomain.coord(:,2));
            zmax = max(obj.meshDomain.coord(:,3));
            xmin = min(obj.meshDomain.coord(:,1));
            ymin = min(obj.meshDomain.coord(:,2));
            zmin = min(obj.meshDomain.coord(:,3));

            % Center and scaling
            a = (xmax - xmin)/2;
            b = (ymax - ymin)/2;
            c = (zmax - zmin)/2;
            x0 = xmin + a;
            y0 = ymin + b;
            z0 = zmin + c;

            % Trilinear shape functions N_i (8 nodes)
            % Local coordinates mapped from global ones
            xi = @(x) (x(1,:,:)-x0)/a;
            eta = @(x) (x(2,:,:)-y0)/b;
            zeta = @(x) (x(3,:,:)-z0)/c;

            % Vector field per node (8 hexahedral nodes)
            N = {
                @(x) 1/8*(1 - xi(x)).*(1 - eta(x)).*(1 - zeta(x));
                @(x) 1/8*(1 + xi(x)).*(1 - eta(x)).*(1 - zeta(x));
                @(x) 1/8*(1 + xi(x)).*(1 + eta(x)).*(1 - zeta(x));
                @(x) 1/8*(1 - xi(x)).*(1 + eta(x)).*(1 - zeta(x));
                @(x) 1/8*(1 - xi(x)).*(1 - eta(x)).*(1 + zeta(x));
                @(x) 1/8*(1 + xi(x)).*(1 - eta(x)).*(1 + zeta(x));
                @(x) 1/8*(1 + xi(x)).*(1 + eta(x)).*(1 + zeta(x));
                @(x) 1/8*(1 - xi(x)).*(1 + eta(x)).*(1 + zeta(x));
                };

            % Create vector fields f{i} = [Nx; Ny; Nz] for each node
            f = cell(1, 8*3);
            for i = 1:8
                f{3*(i-1)+1} = @(x) [N{i}(x); 0*x(1,:,:); 0*x(1,:,:)];  % x-component
                f{3*(i-1)+2} = @(x) [0*x(1,:,:); N{i}(x); 0*x(1,:,:)];  % y-component
                f{3*(i-1)+3} = @(x) [0*x(1,:,:); 0*x(1,:,:); N{i}(x)];  % z-component
            end

            % Create AnalyticalFunction objects
            for i = 1:length(f)
                uD{i} = AnalyticalFunction.create(f{i}, obj.meshDomain.ndim, obj.boundaryMeshJoined);
            end
        end

    end

end
