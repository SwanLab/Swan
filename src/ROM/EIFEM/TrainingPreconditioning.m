classdef TrainingPreconditioning <handle
    
    properties (Access = public)
        mesh
        ndim
        young
        poisson
        meshSource
        bounds
        nelem
        meshFile
        inclusionType
        training
        method
    end

    properties (Access = private)
        material
        levelSet
        geometry
        unfittedMesh
        basisType
        dirac
        coarseOrder
        domainIndices
        nSubdomains

    end
    
    methods (Access = public)
        function obj = TrainingPreconditioning(cParams)
            obj.init(cParams)
        end

        function data = compute(obj)
            obj.createMesh();
            obj.material = obj.createMaterial(obj.mesh);
        
            switch obj.method
                case 'Multiscale'
                    data = obj.computeMultiscale();
                case 'EIFEM'
                    data = obj.computeEIFEM();
            end

        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     methods (Access = private)

         function init(obj,cParams)
            obj.meshSource    = cParams.meshSource;
            obj.ndim          = cParams.ndim;
            obj.young         = cParams.young;
            obj.poisson       = cParams.poisson;
            obj.method        = cParams.method;
            obj.basisType     = cParams.basisType;
            obj.coarseOrder   = cParams.coarseOrder;
            initMesh = sprintf('init%sMesh', obj.meshSource);
            obj.(initMesh)(cParams);

            if strcmp(obj.method, 'EIFEM')
                obj.initEIFEM(cParams);
            end
         end

         function initEIFEM(obj, cParams)
             obj.nSubdomains   = cParams.nSubdomains;
             obj.domainIndices = cParams.domainIndices;
             obj.dirac         = cParams.dirac;
         end

         function initGeneratedMesh(obj, cParams)
             obj.nelem         = cParams.nelem;
             obj.bounds        = cParams.bounds;
             obj.geometry      = cParams.geometry;
             obj.inclusionType = cParams.inclusionType;
         end

         function initImportedMesh(obj, cParams)
             obj.meshFile = cParams.meshFile;
         end

         function createMesh(obj)
            s.meshSource = obj.meshSource;
            s.ndim       = obj.ndim;
            s.nelem      = obj.nelem;
            s.bounds     = obj.bounds;
            s.meshFile   = obj.meshFile;

            m = TrainingMesh(s);
            obj.mesh = m.create();

            if strcmp(obj.meshSource, 'Generated')
                if ~ismember(obj.inclusionType, {'Void', 'Weak'})
                    error('Unknown inclusion type. Use Void or Weak.');
                end

                if strcmp(obj.inclusionType, 'Void')
                    obj.cutMesh();
                end
            end
         end


         function cutMesh(obj)
             bgMesh           = obj.mesh;
             lsFun            = obj.geometry.computeLevelSetFunction(bgMesh);
             obj.levelSet     = lsFun.fValues;
             s.backgroundMesh = bgMesh;
             s.boundaryMesh   = bgMesh.createBoundaryMesh();
             uMesh            = UnfittedMesh(s);
             uMesh.compute(obj.levelSet);
             obj.unfittedMesh = uMesh;
             obj.mesh         = uMesh.createInnerMesh();
         end

         function material = createMaterial(obj,mR)
             % mD = createMeshDomain(nSubdomains,mesh);
             switch obj.meshSource
                 case 'Imported'
                     material = obj.computeElasticProperties(mR);
                 case 'Generated'
                     switch obj.inclusionType
                         case 'Void'
                             material = obj.computeElasticProperties(mR);
                         case 'Weak'
                             s.mesh           = mR;
                             s.geomFun        = obj.geometry;
                             s.bounds         = obj.bounds;
                             s.E              = obj.young;
                             s.nu             = obj.poisson;
                             m = MaterialTraining(s);
                             material = m.create();
                     end
             end
         end

             function m= computeElasticProperties(obj,Mesh)
                 N = obj.mesh.ndim;
                 E = ConstantFunction.create(obj.young,Mesh);
                 nu = ConstantFunction.create(obj.poisson,Mesh);
                 lam = LameParametersConverter;

                 mu     = lam.computeShearFromYoungAndPoisson(E, nu);  mu = Expand(mu,4);
                 lambda = lam.computeLambdaFromYoungAndPoisson(E, nu,N);  lambda = Expand(lambda,4);
                 I      = ConstantFunction.create(eye4D(N),Mesh);
                 IxI    = ConstantFunction.create(kronEye(N),Mesh);

                 m= 2*mu.*I + lambda.*IxI;
             end


         function data=computeMultiscale(obj)
             switch obj.basisType
                 case 'Continuous'
                     s.bMesh        = obj.mesh.createSingleBoundaryMesh();
                     s.lambdaFun  = LagrangianFunction.create(s.bMesh,obj.mesh.ndim, 'P1');
                 case 'Discontinuous'
                     s.bMesh        = obj.mesh.createBoundaryMesh();
                     for i=1:numel(s.bMesh)
                         s.lambdaFun{i}=LagrangianFunction.create(s.bMesh{i}.mesh,obj.mesh.ndim, 'P1');
                     end
             end
             s.mesh          = obj.mesh;
             s.material      = obj.material;
             s.type          = obj.basisType;
             s.uFun          = LagrangianFunction.create(s.mesh, s.mesh.ndim, 'P1');
             s.dirichletFun  = obj.createDirichletFunction(s.bMesh);
             e  = ElasticHarmonicExtension(s);
             [T,~,~,Kcoarse] = e.solve();
             data.T       = T;
             data.Kcoarse = Kcoarse;
             data.mesh    = obj.mesh;

         end

         function data = computeEIFEM(obj)
             s.mesh           = obj.mesh;
             s.domainIndices  = obj.domainIndices;
             s.nSubdomains    = obj.nSubdomains;
             s.type           = obj.basisType;
             s.Coarseorder    = obj.coarseOrder;
             e= EIFEMTraining(s);
             matDomain = obj.createMaterial(e.meshDomain);
             Eif              = e.train(matDomain);
             Eif.dirac        = obj.dirac;
             Eif.type         = obj.basisType;
             Eif.Coarseorder  = 1;
             Eif.material     = obj.material;

             z = OfflineDataProcessor(Eif);
             EIFEoper      = z.computeROMbasis();
             data.T        = EIFEoper.T;
             data.mesh     = obj.mesh;
             data.Kcoarse  = EIFEoper.Kcoarse;
         end


         function dF = createDirichletFunction(obj,bMesh)
             s.mesh = bMesh;
             s.order= obj.coarseOrder;
             s.type = obj.basisType;
             cf = CoarseFunctions(s);
             dF = cf.getAnalytical();
         end



     end
end

