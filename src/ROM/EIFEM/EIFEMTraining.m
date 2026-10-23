classdef EIFEMTraining < handle

    properties (GetAccess = public, SetAccess = private)
        refMesh
        Coarseorder
    end

    properties (Access = private)
        meshDomain
        DDdofManager
        domainIndices
        tolSameNode
        cellMeshes
        nSubdomains
        geometryType
        levelSet
        unfittedMesh
        type
        young   
        poisson  
    end

    methods (Access = public)

        function obj = EIFEMTraining(cParams)
            obj.init(cParams);
        end

        function data = train(obj)
            obj.createTrainingMesh();
            bMesh     = obj.createBoundaryMesh();
            material  = obj.createMaterial();

            uFun      = obj.createDisplacementFun();
            lambdaFun = obj.createLagrangeMultiplier(bMesh);
            coarseFun = obj.createCoarseFunctions(bMesh);

            [u, ~, K, ~] = obj.solveHarmonicExtension(obj.meshDomain,bMesh,material,uFun,lambdaFun,coarseFun,obj.type);
            
            [data.uSbd, data.LHSsbd] = obj.extractDomainData(u, K);
            data.mesh = obj.refMesh;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.refMesh       = cParams.mesh;
            obj.nSubdomains   = cParams.nSubdomains;
            obj.domainIndices = cParams.domainIndices;
            obj.tolSameNode   = 1e-11;
            obj.Coarseorder   = cParams.Coarseorder;
            obj.type          = cParams.type;
            obj.young         = cParams.young;
            obj.poisson       = cParams.poisson;
        end

        function material = createMaterial(obj)
            obj.young   = 1;
            obj.poisson = 1/3;  

            N = obj.meshDomain.ndim;
            E = ConstantFunction.create(obj.young,obj.meshDomain);
            nu = ConstantFunction.create(obj.poisson,obj.meshDomain); 
            lam = LameParametersConverter;

            mu     = lam.computeShearFromYoungAndPoisson(E,nu);  mu = Expand(mu,4);
            lambda = lam.computeLambdaFromYoungAndPoisson(E,nu,N);  lambda = Expand(lambda,4);
            I      = ConstantFunction.create(eye4D(N),obj.meshDomain);
            IxI    = ConstantFunction.create(kronEye(N),obj.meshDomain);
            
            material = 2*mu.*I + lambda.*IxI;
        
        end

        function bMesh = createBoundaryMesh(obj)
            switch obj.type
                case 'Continuous'
                    bMesh = obj.meshDomain.createSingleBoundaryMesh();
                case 'Discontinuous'
                    bMesh = obj.meshDomain.createBoundaryMesh();
            end
        end

        function lambdaFun = createLagrangeMultiplier(obj,bMesh)
            switch obj.type
                case 'Continuous'
                    lambdaFun = LagrangianFunction.create(bMesh, obj.meshDomain.ndim, 'P1');
                case 'Discontinuous'
                    lambdaFun = cell(1, numel(bMesh));
                    for i = 1:numel(bMesh)
                        lambdaFun{i} = LagrangianFunction.create(bMesh{i}.mesh, obj.meshDomain.ndim, 'P1');
                    end
            end
        end

        function uFun = createDisplacementFun(obj)
             uFun      = LagrangianFunction.create(obj.meshDomain, obj.meshDomain.ndim, 'P1');
        end

        function createTrainingMesh(obj)
            if sum(obj.nSubdomains > 1) >= 1
                bS = obj.refMesh.createBoundaryMesh();
                [mD, mSb, iC, lG, iCR, ~] = obj.createMeshDomain(obj.refMesh);
                obj.cellMeshes   = mSb;
                obj.meshDomain   = mD;
                obj.DDdofManager = obj.createDomainDecompositionDofManager(iC, lG, bS, obj.refMesh, iCR);
            else
                obj.cellMeshes = {obj.refMesh};
                obj.meshDomain = obj.refMesh;
            end
        end

        function [mD, mSb, iC, lG, iCR, discMesh] = createMeshDomain(obj, mR)
            s.nsubdomains   = obj.nSubdomains;
            s.meshReference = mR;
            s.tolSameNode   = obj.tolSameNode;
            s.geometryType  = obj.geometryType;
            m               = MeshCreatorFromRVE.create(s);
            [mD, mSb, iC, ~, lG, iCR, discMesh] = m.create();
        end

        function [u, lhs] = extractDomainData(obj, uC, LHS)
            if sum(obj.nSubdomains > 1) >= 1
                u   = obj.extractDomainDisplacements(uC);
                lhs = 1;
            else
                u   = full(uC);
                lhs = LHS;
            end
        end

        function u = extractDomainDisplacements(obj, uC)
            ntest = size(uC, 2);
            ind   = obj.getFlatSubdomainIndex();
            u     = zeros(obj.refMesh.nnodes*obj.refMesh.ndim, ntest); 
            
            for i = 1:ntest
                uD = obj.DDdofManager.global2local(uC(:, i));
                u(:, i) = uD(:, ind);
            end
        end

        function lhs = extractDomainLHS(obj, LHS)
            lhs = obj.DDdofManager.global2localMatrix(LHS);
            ind = obj.getFlatSubdomainIndex();
            lhs = lhs(:, :, ind);
        end

        function ind = getFlatSubdomainIndex(obj)
            % Encapsulates 2D/3D subdomain index mapping
            if obj.refMesh.ndim == 2
                ind = (obj.domainIndices(1) - 1) * obj.nSubdomains(1) + obj.domainIndices(2);
            else
                ind = obj.nSubdomains(1) * obj.nSubdomains(2) * (obj.domainIndices(3) - 1) ...
                    + obj.nSubdomains(1) * (obj.domainIndices(2) - 1) + obj.domainIndices(1);
            end
        end

        function d = createDomainDecompositionDofManager(obj, iC, lG, bS, mR, iCR)
            s.nSubdomains            = obj.nSubdomains;
            s.interfaceConnec        = iC;
            s.interfaceConnecReshaped = iCR;
            s.locGlobConnec          = lG;
            s.nBoundaryNodes         = bS{1}.mesh.nnodes;
            s.nReferenceNodes        = mR.nnodes;
            s.nNodes                 = obj.meshDomain.nnodes;
            s.nDimf                  = obj.meshDomain.ndim;

            % Dynamic instantiation based on dimension
            dimSuffix    = ternary(obj.refMesh.ndim == 3, '3D', '');
            managerClass = ['DomainDecompositionDofManager' dimSuffix];
            d            = feval(managerClass, s);
        end

        function dF = createCoarseFunctions(obj, bMesh)
            s.mesh  = bMesh;
            s.order = obj.Coarseorder;
            s.type  = obj.type;
            cf      = CoarseFunctions(s);
            dF      = cf.getAnalytical();
        end

        function [u, lambda, K, Kc] = solveHarmonicExtension(obj,mesh,bMesh,material,uFun,lambdaFun,coarseFun,type)
            s.mesh         = mesh;
            s.bMesh        = bMesh;
            s.lambdaFun    = lambdaFun;
            s.uFun         = uFun;
            s.material     = material;
            s.dirichletFun = coarseFun;
            s.type         = type;

            elasExt            = ElasticHarmonicExtension(s);
            [u, lambda, K, Kc] = elasExt.solve();
        end

    end

end

function val = ternary(cond, trueVal, falseVal)
    if cond, val = trueVal; else, val = falseVal; end
end