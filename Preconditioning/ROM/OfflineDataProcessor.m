classdef OfflineDataProcessor < handle

    properties (GetAccess = public , SetAccess = private)

    end
    properties (Access = private)
        mesh
        boundaryMeshJoined
        localGlobalConnecBd
        LHS
        Coarseorder
        fValuesTraining
        RigidBodyFun
        DeformationalFun
        material

        fileNameData

    end

    methods (Access = public)

        function obj = OfflineDataProcessor(cParams)
            obj.init(cParams)
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

            [Vfun,Vfun2] = obj.createInterfaceModesFun(bMesh);


            uDefFunBd  = obj.restrictToBoundary(uDefFun,bMesh);
            RBFunBd    = obj.restrictToBoundary(uRBfun(1),bMesh); %only the first bevause we just want the basis!
            LMDefFunBd = obj.restrictToBoundary(LMDefFun,bMesh);

            Adr = obj.computeBoundaryModalMassMatrix(uDefFunBd,RBFunBd);
            Arr = obj.computeBoundaryModalMassMatrix(RBFunBd,RBFunBd);
            Add = obj.computeBoundaryModalMassMatrix(uDefFunBd,LMDefFunBd);
            Ldv = obj.computeBoundaryModalMassMatrix(LMDefFunBd,Vfun);
            Lrv = obj.computeBoundaryModalMassMatrix(RBFunBd,Vfun);            

            RBFun = uRBfun(1);
            RBFun2.basisFunctions{1} = project(RBFun.basisFunctions{1},'P1');
            RBFun2.basisFunctions{2} = project(RBFun.basisFunctions{2},'P1');
            RBFun2.basisFunctions{3} = project(RBFun.basisFunctions{3},'P1');

            uDefFunBd2 = obj.restrictToBoundary2(uDefFun);            
            LMDefFunBd2 = obj.restrictToBoundary2(LMDefFun);
            RBFunBd2    = obj.restrictToBoundary2(RBFun2);

            
            Adr2 = obj.computeBoundaryModalMassMatrix2(uDefFunBd2,RBFunBd2);           
            Arr2 = obj.computeBoundaryModalMassMatrix2(RBFunBd2,RBFunBd2);
            Add2 = obj.computeBoundaryModalMassMatrix2(uDefFunBd2,LMDefFunBd2);          
            Ldv2 = obj.computeBoundaryModalMassMatrix2(LMDefFunBd2,Vfun2);
            Lrv2 = obj.computeBoundaryModalMassMatrix2(RBFunBd2,Vfun2);            



            Ud = PhiD*(Add'\Ldv);
            Ur = PhiR*inv(Arr')*(Lrv - Adr'*(Add'\Ldv));
            U  = Ur+Ud;


            nB = 8;
            Kdd = PhiD'*obj.LHS*PhiD;
            %
            nld = LMDefFun.nbasis;
            nlr = uRBfun(1).nbasis;            
            Zrd = zeros(nlr,nld);
            Zd = zeros(nld,nB);
            Zr = zeros(nlr,nB);
            Zrr = zeros(nlr,nlr);

            Keif = [Kdd Zrd';Zrd Zrr];
            C    = [Adr Add;...
                    Arr Zrd];
 
            

            Z  = zeros(nld+nlr,nld+nlr);
%
            LHS = [Keif C; C.' Z];
            Lug = [Lrv;Ldv]*eye(8);
            RHS = [Zd;Zr;Lug];
%
            x = LHS\RHS;
            %
            
            uEifD = x(1:nld,:);
            uEifR = x(nld+1:nld+nlr,:);
            
            U2 = PhiD*uEifD + PhiR*uEifR;



            % s.mesh         = obj.mesh;
            % s.uFun         = cParams.uFun;
            % s.lambdaFun    = cParams.lambdaFun;
            % s.material     = cParams.material;
            % s.dirichletFun = cParams.dirichletFun;
            % s.localGlobalConnecBd = cParams.localGlobalConnecBd;
            % 
            % 
            % 
            % %ElasticHarmonicExtension
            % 

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

        function init(obj,cParams)
            obj.mesh            = cParams.mesh;
            obj.fValuesTraining = cParams.uSbd;
            obj.LHS             = cParams.LHSsbd;
            obj.Coarseorder     = cParams.Coarseorder;
            obj.material        = cParams.material;
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

        function BdFun = restrictToBoundary2(obj,fun)
            nboundary = numel(fun.basisFunctions);
            for i = 1:nboundary
                BdFun{i} = fun.basisFunctions{i}.restrictToBoundary();
            end
        end        

        function BdFun = restrictToBoundary(obj,fun,bMesh)
            nboundary = size(bMesh,1);
            for i = 1:nboundary
                mesh  = bMesh{i};
                BdFun{i} = fun.restrictBasisToBoundaryMesh(mesh);
            end
        end

        function M = computeBoundaryModalMassMatrix2(obj,test,trial)
            nTest  = numel(test);
            nTrial = numel(trial);
            M = zeros(nTest,nTrial);
            quadOrder = 2;
            for i = 1:nTest
                for j = 1:nTrial
                    v = test{i};
                    u = trial{j};
                    int = DP(u,v);
                    M(i,j) = Integrator.compute(int,int.mesh,quadOrder);
                end
            end
        end

        function M = computeBoundaryModalMassMatrix(obj,test,trial)
            nTest  = test{1}.nbasis;
            nTrial = trial{1}.nbasis;
            M = zeros(nTest,nTrial);
            nbd = size(test,2);
            quadOrder = 2;
            for ibd = 1:nbd
                for i = 1:nTest
                    for j = 1:nTrial
                        v = test{ibd}.basisFunctions{i};
                        u = trial{ibd}.basisFunctions{j};
                        int = DP(u,v);
                        M(i,j) = M(i,j) + Integrator.compute(int,int.mesh,quadOrder);
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

        function [Vfun,Vfun2] = createInterfaceModesFun(obj,bMesh)
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

            boundMesh=obj.mesh.createSingleBoundaryMesh();
                  nfun = length(f);
              for k = 1:nfun
                  uD{k} = AnalyticalFunction.create(f{k}, boundMesh);
              end
            f = uD;
            bMesh= bMesh([1 3 2 4]); %Reorder to get anticlock wise boundary functions
            s.mesh=boundMesh;
            s.type='quad';
            cf=CoarseFunctions(s);
            f2=cf.compute();
            f3=InterfaceFunctions(boundMesh);

            %CoarseFunction 
            %Netejar Coarse funciton

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

            bMesh2   = obj.mesh.createSingleBoundaryMesh();            
            Vfun2=cell(1,8);
            for i=1:nfun
                Vfun2{i}  = AnalyticalFunction.create(f{i},bMesh2);
            end            
        end


        function [LHS,u] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
            LHS = obj.computeLHS(u);
%             RHS = obj.computeRHS(u,dLambda);
        end

        function K  = computeLHS(obj,u)          
            C = obj.material;
            K = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),u,u,obj.mesh,'Domain',2);
        end

        function M  = computeM(obj,u)          
            M = IntegrateLHS(@(u,v) DD(v,u),u,u,obj.mesh,'Domain',2);
        end


    end

end


