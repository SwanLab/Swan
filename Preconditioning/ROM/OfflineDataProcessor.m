classdef OfflineDataProcessor < handle

    properties (GetAccess = public , SetAccess = private)

    end
    properties (Access = private)
        mesh
        boundaryMeshJoined
        localGlobalConnecBd
        LHS
        E
        nu
        Coarseorder
        fValuesTraining
        RigidBodyFun
        DeformationalFun

        fileNameData

    end

    methods (Access = public)

        function obj = OfflineDataProcessor(data)
            obj.init(data)
        end

        function EIFEoper = computeROMbasis(obj)
            obj.LHS = createElasticProblem(obj);

            uFun         = obj.createDispFun();
            uRBfun       = obj.projectToRigidBodyFun(uFun);
            uDEFSpaceFun = obj.projectToDeformationalSpace(uFun,uRBfun);
            PhiD         = obj.computeDeformationalModes(uDEFSpaceFun);
            uDefFun      = obj.createDeformationalFunction(PhiD);

            PsiD       = obj.computeSelfEquilibratedLagrange(PhiD);
            LMDefFun   = obj.createDeformationalFunction(PsiD);

            PhiR       = obj.getRigidBodyModes(uRBfun(1));

%             u = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
%             M = IntegrateLHS(@(u,v) DP(v,u),u,u,obj.mesh,'Domain',2);
%             defSnap = obj.fValuesTraining - PhiR*inv(PhiR'*M*PhiR)*(PhiR'*M*obj.fValuesTraining);

            bMesh = obj.mesh.createBoundaryMesh();

            Vfun = obj.createInterfaceModesFun(bMesh);

            uDefFunBd  = obj.restrictToBoundary(uDefFun,bMesh);
            RBFunBd    = obj.restrictToBoundary(uRBfun(1),bMesh); %only the first bevause we just want the basis!
            LMDefFunBd = obj.restrictToBoundary(LMDefFun,bMesh);

            Adr = obj.computeBoundaryModalMassMatrix(uDefFunBd,RBFunBd);
            Arr = obj.computeBoundaryModalMassMatrix(RBFunBd,RBFunBd);
            %             Add = obj.computeBoundaryModalMassMatrix(uDefFunBd,LMDefFunBd);
            Add = obj.computeBoundaryModalMassMatrixDirac(uDefFunBd,LMDefFunBd);
            %             Add = PhiD'*PsiD;
            % Add with dirac integrator and PhiD'*PsiD are the same!
%             Ldv = obj.computeBoundaryModalMassMatrix(LMDefFunBd,Vfun);
            Ldv = obj.computeBoundaryModalMassMatrixDirac(LMDefFunBd,Vfun);
            Lrv = obj.computeBoundaryModalMassMatrix(RBFunBd,Vfun);

            Ud = PhiD*(Add'\Ldv);
            Ur = PhiR*inv(Arr')*(Lrv - Adr'*(Add'\Ldv));
            U  = Ur+Ud;


            %ndof = 8;
            %Kdd = PhiD'*obj.LHS*phiD;
            %
            %Keif = [Kdd Zdr;Zrd Zrr];
            %C    = [Adr Add;...
            %        Arr S];
%
            %LHS = [K C; C.' Z];
            %Lug = [Lrv;Ldv]*eye(8);
            %RHS = -[Zd;Zr;Lug];
%
            %x = LHS\RHS;
            %
            %U2 = x(1:ndof*2,1);


            Kcoarse = Ud'*obj.LHS*Ud;

            EIFEoper.Kcoarse = Kcoarse;
            EIFEoper.U       = U;
            EIFEoper.Urb     = Ur;
            EIFEoper.Udef    = Ud;
            EIFEoper.PhiD    = PhiD;
            EIFEoper.PhiR    = PhiR;
            EIFEoper.Kfine   = obj.LHS;
        end

    end

    methods (Access = private)

        function init(obj,data)
            obj.mesh            = data.mesh;
            obj.fValuesTraining = data.uSbd;
            obj.LHS             = data.LHSsbd;
            obj.E               = data.E;
            obj.nu              = data.nu;
            obj.Coarseorder     = data.Coarseorder;
        end

        function uFun = createDispFun(obj)
            ntest  = size(obj.fValuesTraining,2);
            s.mesh  = obj.mesh;
            s.ndimf = 2;
            s.order = 'P1';
            for i = 1:ntest
                fValues   = obj.fValuesTraining(:,i);
                s.fValues = reshape(fValues,2,[])' ;
                uFun(i)   = LagrangianFunction(s);
            end
        end

        function RBfun = projectToRigidBodyFun(obj,uFun)
            refPoint = (max(obj.mesh.coord)+min(obj.mesh.coord))/2;
            ntest    = size(uFun,2);
            for i = 1:ntest
                RBfun(i) = project(uFun(i),'RigidBody',refPoint);
            end
        end

        function DEFfun = projectToDeformationalSpace(obj,uFun,RBfun)
            ntest    = size(uFun,2);
            for i = 1:ntest
                Df = uFun(i)-RBfun(i);
%                 s.operation = @(xV) uFun(i).evaluate(xV) - RBfun(i).evaluate(xV);
%                 s.ndimf   = uFun.ndimf;
%                 s.mesh    = uFun.mesh;
%                 Df        = DomainFunction(s);
                DEFfun(i) = project(Df,'P1');
            end
        end

        function phiD = computeDeformationalModes(obj,DEFfun)
            val = obj.constructSnapshotMatrix(DEFfun);
            [U,S,~] = svd(val,"econ");
            tol = 1e-6;
            r = sum(diag(S) > tol);  % Count significant singular values
            phiD = U(:, 1:r);% Keep only the first r singular values
        end

        function fV = constructSnapshotMatrix(obj,fun)
            ntest    = size(fun,2);
            fV =[];
            for i = 1:ntest
                val = fun(i).fValues;
                val = reshape(val',[],1);
                fV  = [fV val];
            end
        end

        function DefFun = createDeformationalFunction(obj,phiCoeff)
            phiCoeff = num2cell(phiCoeff,1);
            phiCoeff = cellfun(@(x) reshape(x', 2,[])', phiCoeff, 'UniformOutput', false);
            functionType = {'P1'};
            functionType = repelem(functionType,size(phiCoeff,2));

            DefFun = ModalFunction.create(obj.mesh,phiCoeff,functionType);
        end

        function Psid = computeSelfEquilibratedLagrange(obj,phid)
            Psid = obj.LHS*phid;
        end

        function PhiR = getRigidBodyModes(obj,RBfun)
            nbasis = RBfun.nbasis;
            for i = 1:nbasis
                fun = project(RBfun.basisFunctions{i},'P1');
                fValues = fun.fValues;
                fValues = reshape(fValues',[],1);
                PhiR(:,i) = fValues;
            end

        end

        function BdFun = restrictToBoundary(obj,fun,bMesh)
            nboundary = size(bMesh,1);
            for i = 1:nboundary
                mesh  = bMesh{i};
                BdFun{i} = fun.restrictBasisToBoundaryMesh(mesh);
            end
        end

        function M = computeBoundaryModalMassMatrix(obj,test,trial)
            nbasistest  = test{1}.nbasis;
            nbasistrial = trial{1}.nbasis;
            M   = zeros(nbasistest,nbasistrial);
            nFlds = test{1}.ndimf;
            nbd = size(test,2);
            for ibd = 1:nbd
                bMesh = test{ibd}.mesh;
                quad  = obj.createQuadrature(bMesh);
                xV    = quad.posgp;
                dV    = bMesh.computeDvolume(quad);
                basisTest   = test{ibd}.evaluateBasisFunctions(xV);
                basisTrial  = trial{ibd}.evaluateBasisFunctions(xV);
                for i = 1:nbasistest
                    for j = 1:nbasistrial
                        for iField = 1:nFlds
                            basisProd = squeeze(basisTest{i}(iField,:,:).*basisTrial{j}(iField,:,:));
                            M(i,j)  = M(i,j) + sum(basisProd.*dV,'all');
                        end
                    end
                end
            end
        end

        function M = computeBoundaryModalMassMatrixDirac(obj,test,trial)
            nbasistest  = test{1}.nbasis;
            nbasistrial = trial{1}.nbasis;
            M   = zeros(nbasistest,nbasistrial);
            nFlds = test{1}.ndimf;
            nbd = size(test,2);
            for ibd = 1:nbd
                bMesh = test{ibd}.mesh;
                for i = 1:nbasistest
                    basisTest   = test{ibd}.basisFunctions{i}.fValues;
                    for j = 1:nbasistrial
                        basisTrial  = trial{ibd}.basisFunctions{j}.fValues;
                        for iField = 1:nFlds
                            basisProd = squeeze(basisTest(:,iField).*basisTrial(:,iField));
                            M(i,j)  = M(i,j) + sum(basisProd,'all');
                        end
                    end
                end
            end
        end

        function q = createQuadrature(obj,mesh)
            order = 2;
            q = Quadrature.create(mesh,order);
        end

        function Vfun = createInterfaceModesFun(obj,bMesh)
            xmax = max(obj.mesh.coord(:,1));
            ymax = max(obj.mesh.coord(:,2));
            xmin = min(obj.mesh.coord(:,1));
            ymin = min(obj.mesh.coord(:,2));
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

            f     = {f1x f1y  f2x f2y f3x f3y f4x f4y}; %
            nfun = size(f,2);
            nbd = size(bMesh,1);
            for ibd=1:nbd
                mesh = bMesh{ibd}.mesh;
                for i=1:nfun
                    uD        = AnalyticalFunction.create(f{i},mesh);
                    uD        = project(uD,'P1');
                    VCoeff{i} = uD.fValues;
                end
                %                 VCoeff = cellfun(@(x) reshape(x', 2,[])', VCoeff, 'UniformOutput', false);
                functionType = {'P1'};
                functionType = repelem(functionType,size(VCoeff,2));
                Vfun{ibd} = ModalFunction.create(mesh,VCoeff,functionType);
            end
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
            E1  = 1;
            E2 = E1/1000;
            nu = 1/3;
%            young   = ConstantFunction.create(obj.E,mesh);
%            poisson = ConstantFunction.create(obj.nu,mesh);
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

        function [LHS,u] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
            LHS = obj.computeLHS(u);
%             RHS = obj.computeRHS(u,dLambda);
        end

        function K  = computeLHS(obj,u)          
            C = obj.createMaterial(obj.mesh);
            K = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),u,u,obj.mesh,'Domain',2);
        end

        function M  = computeM(obj,u)          
            M = IntegrateLHS(@(u,v) DD(v,u),u,u,obj.mesh,'Domain',2);
        end


    end

end


