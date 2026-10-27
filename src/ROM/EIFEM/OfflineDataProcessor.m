classdef OfflineDataProcessor < handle

    properties (GetAccess = public , SetAccess = private)

    end
    properties (Access = private)
        mesh
        LHS
        Coarseorder
        fValuesTraining
        young   
        poisson  
        dirac
        type
    end

    methods (Access = public)

        function obj = OfflineDataProcessor(cParams)
            obj.init(cParams)
        end

        function EIFEoper = computeROMbasis(obj)
            obj.LHS              = obj.computeLHS();
            bMesh                = obj.mesh.createBoundaryMesh();
            
            [uDef,uRB,PhiD,PhiR] = obj.createReducedDispBasis();
            [LMdef,LMrb]         = obj.createLagrangeMultBasis(PhiD,uRB,bMesh);
            Vfun                 = obj.createInterfaceModesFun(bMesh);

            [Add,Adr,Arr,Ldv,Lrv] = obj.computeReducedMatrices(uDef,uRB,LMdef,LMrb,Vfun,bMesh);           
            [T,Td,Tr]             = obj.computeInterscaleOperator(PhiD,PhiR,Add,Adr,Arr,Ldv,Lrv);


            % uFun         = obj.createDispFun();
            % uRBfun       = obj.projectToRigidBodyFun(uFun);
            % uDEFSpaceFun = obj.projectToDeformationalSpace(uFun,uRBfun);
            % PhiD         = obj.computeDeformationalModes(uDEFSpaceFun);
            % uDefFun      = obj.createDeformationalFunction(PhiD);
          
            % PsiD       = obj.computeSelfEquilibratedLagrange(PhiD);
            % LMDefFun   = obj.createDeformationalFunction(PsiD);
            % PhiR       = obj.getRigidBodyModes(uRBfun(1));       
            % Vfun       = obj.createInterfaceModesFun(bMesh);

            % uDefFunBd  = obj.restrictToBoundary(uDef,bMesh);
            % uRBFunBd   = obj.restrictToBoundary(uRBfun(1),bMesh); %only the first bevause we just want the basis!
            % % LMDefFunBd = obj.restrictToBoundary(LMDefFun,bMesh);
            % 
            % Adr = obj.computeBoundaryModalMassMatrix(uDefFunBd,LMrb);
            % Arr = obj.computeBoundaryModalMassMatrix(uRBFunBd,LMrb);
            % Lrv = obj.computeBoundaryModalMassMatrix(LMrb,Vfun); 
            % if obj.dirac==true
            %     Add = obj.computeBoundaryModalMassMatrixDirac(uDefFunBd,LMDefFunBd);
            %     Ldv = obj.computeBoundaryModalMassMatrixDirac(LMdef,Vfun);
            % else
            %     Add = obj.computeBoundaryModalMassMatrix(uDefFunBd,LMDefFunBd);
            %     Ldv = obj.computeBoundaryModalMassMatrix(LMdef,Vfun);
            % end

            
            Kcoarse = T'*obj.LHS*T;

            EIFEoper.Kcoarse = Kcoarse;
            EIFEoper.T       = T;
            EIFEoper.Trb     = Tr;
            EIFEoper.Tdef    = Td;
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
            obj.young           = cParams.young;
            obj.poisson         = cParams.poisson; 
            obj.type            = cParams.type;

            if isfield(cParams,'dirac')
                obj.dirac = cParams.dirac;
            else
                obj.dirac = true;
            end

        end

        function uFun = createDispFun(obj)
            ntest  = size(obj.fValuesTraining,2);
            s.mesh  = obj.mesh;
            s.ndimf = obj.mesh.ndim;
            s.order = 'P1';
            for i = 1:ntest
                fValues   = obj.fValuesTraining(:,i);
                s.fValues = reshape(fValues,obj.mesh.ndim,[])' ;
                uFun(i)   = LagrangianFunction(s);
            end
        end

        function RBfun = projectToRigidBodyFun(obj,uFun)
            refPoint = (max(obj.mesh.coord)+min(obj.mesh.coord))/2;
            ntest    = size(uFun,2);
            for i = 1:ntest
                RBfun(i) = uFun(i).project('RigidBody',refPoint);
                % RBfun(i) = project(uFun(i),);
            end
        end

        function DEFfun = projectToDeformationalSpace(~,uFun,RBfun)
            ntest = size(uFun,2);
            for i = 1:ntest
                Df = uFun(i)-RBfun(i);
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

        function fV = constructSnapshotMatrix(~,fun)
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
            phiCoeff = cellfun(@(x) reshape(x', obj.mesh.ndim ,[])', phiCoeff, 'UniformOutput', false);
            functionType = {'P1'};
            functionType = repelem(functionType,size(phiCoeff,2));

            DefFun = ModalFunction.create(obj.mesh,phiCoeff,functionType);
        end

        function Psid = computeSelfEquilibratedLagrange(obj,phid)
            Psid = obj.LHS*phid;
        end

        function PhiR = getRigidBodyModes(~,RBfun)
            nbasis = RBfun.nbasis;
            for i = 1:nbasis
                fun = project(RBfun.basisFunctions{i},'P1');
                fValues = fun.fValues;
                fValues = reshape(fValues',[],1);
                PhiR(:,i) = fValues;
            end
        end

        function BdFun = restrictToBoundary2(~,fun)
            nboundary = numel(fun.basisFunctions);
            for i = 1:nboundary
                BdFun{i} = fun.basisFunctions{i}.restrictToBoundary();
            end
        end        

        function BdFun = restrictToBoundary(~,fun,bMesh)
            nboundary = size(bMesh,1);
            for i = 1:nboundary
                Mesh  = bMesh{i};
                BdFun{i} = fun.restrictBasisToBoundaryMesh(Mesh);
            end
        end

        function M = computeBoundaryModalMassMatrix2(~,test,trial)
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

        function M = computeBoundaryModalMassMatrix(~,test,trial)
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

        function M = computeBoundaryModalMassMatrixDirac(~,test,trial)
            nbasistest  = test{1}.nbasis;
            nbasistrial = trial{1}.nbasis;
            M   = zeros(nbasistest,nbasistrial);
            nFlds = test{1}.ndimf;
            nbd = size(test,2);
            for ibd = 1:nbd
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


        function Vfun = createInterfaceModesFun(obj,bMesh)
            switch obj.type
                case 'Continuous'
                    joinedBMesh= obj.mesh.createSingleBoundaryMesh();
                    s.mesh     = joinedBMesh;
                    s.type     = obj.type;
                    s.order    = obj.Coarseorder;
                    cf         = CoarseFunctions(s);
                    f          = cf.compute();
                    nfun = size(f,2);
                    nbd  = size(bMesh,1);
                    Vfun = cell(1,nbd);
                    for ibd = 1:nbd
                        Mesh = bMesh{ibd}.mesh;
                        VCoeff = cell(1, nfun);
                        ftype = cell(1, nfun);
                        for i = 1:nfun
                            uD = AnalyticalFunction.create(f{i}, Mesh);
                            uD = project(uD, 'P1');
                            VCoeff{i} = uD.fValues;
                            ftype{i} = 'P1';
                        end
                        Vfun{ibd} = ModalFunction.create(Mesh, VCoeff, ftype);
                    end
           
                case 'Discontinuous'
                    s.mesh  = bMesh;
                    s.type  = obj.type;
                    s.order = obj.Coarseorder;
                    cf      = CoarseFunctions(s);
                    f       = cf.compute();
                    nfun    = size(f,2);
                    nbd     = size(bMesh,1);
                    nPerSeg = nfun/nbd;
                    Vfun = cell(1,nbd);
                    k=1;
                    for ibd = 1:nbd
                        Mesh = bMesh{ibd}.mesh;
                        VCoeff = cell(1, nfun);
                        ftype = cell(1, nfun);
                        for i = 1:nfun
                            VCoeff{i} = zeros(Mesh.nnodes,Mesh.ndim);
                            ftype{i}  = 'P1';
                        end
                        for i = 1:nPerSeg
                            uD = AnalyticalFunction.create(f{k}, Mesh);
                            uD = project(uD, 'P1');
                            VCoeff{k} = uD.fValues;
                            ftype{k} = 'P1';
                            k=k+1;
                        end
                        Vfun{ibd} = ModalFunction.create(Mesh, VCoeff, ftype);
                    end
            end
        end

        function [uDef,uRB,PhiD,PhiR] = createReducedDispBasis(obj)
            uFun         = obj.createDispFun();
            uRB          = obj.projectToRigidBodyFun(uFun);
            uDEFSpaceFun = obj.projectToDeformationalSpace(uFun,uRB);
            PhiD         = obj.computeDeformationalModes(uDEFSpaceFun);
            uDef         = obj.createDeformationalFunction(PhiD);
            PhiR         = obj.getRigidBodyModes(uRB(1));    
        end

        function [LMdef,LMrb] = createLagrangeMultBasis(obj,PhiD,uRB,bMesh)
            PsiD       = obj.computeSelfEquilibratedLagrange(PhiD);
            LMdef      = obj.createDeformationalFunction(PsiD);
            LMdef      = obj.restrictToBoundary(LMdef,bMesh);
            LMrb       = obj.restrictToBoundary(uRB(1),bMesh);
        end

        function [Add,Adr,Arr,Ldv,Lrv] = computeReducedMatrices(obj,uDef,uRB,LMdef,LMrb,Vfun,bMesh)
            uDefFunBd  = obj.restrictToBoundary(uDef,bMesh);
            uRBFunBd   = obj.restrictToBoundary(uRB(1),bMesh); %only the first bevause we just want the basis!
            % LMDefFunBd = obj.restrictToBoundary(LMDefFun,bMesh);

            Adr = obj.computeBoundaryModalMassMatrix(uDefFunBd,LMrb);
            Arr = obj.computeBoundaryModalMassMatrix(uRBFunBd,LMrb);
            Lrv = obj.computeBoundaryModalMassMatrix(LMrb,Vfun);
            if obj.dirac==true
                Add = obj.computeBoundaryModalMassMatrixDirac(uDefFunBd,LMdef);
                Ldv = obj.computeBoundaryModalMassMatrixDirac(LMdef,Vfun);
            else
                Add = obj.computeBoundaryModalMassMatrix(uDefFunBd,LMdef);
                Ldv = obj.computeBoundaryModalMassMatrix(LMdef,Vfun);
            end
        end

        function K  = computeLHS(obj)          
            C = obj.createMaterial();
            u = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
            K = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),u,u,obj.mesh,'Domain',2);
        end

        function material = createMaterial(obj)
            obj.young   = 1;
            obj.poisson = 1/3;  

            N = obj.mesh.ndim;
            E = ConstantFunction.create(obj.young,obj.mesh);
            nu = ConstantFunction.create(obj.poisson,obj.mesh); 
            lam = LameParametersConverter;

            mu     = lam.computeShearFromYoungAndPoisson(E,nu);  mu = Expand(mu,4);
            lambda = lam.computeLambdaFromYoungAndPoisson(E,nu,N);  lambda = Expand(lambda,4);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            
            material = 2*mu.*I + lambda.*IxI;
        
        end

        function M  = computeM(obj,u)          
            M = IntegrateLHS(@(u,v) DD(v,u),u,u,obj.mesh,'Domain',2);
        end

        function [T,Td,Tr] = computeInterscaleOperator(~,PhiD,PhiR,Add,Adr,Arr,Ldv,Lrv)
            Td = PhiD*(Add'\Ldv);
            Tr = PhiR*inv(Arr')*(Lrv - Adr'*(Add'\Ldv));
            T  = Tr+Td;
        end

    end

end


