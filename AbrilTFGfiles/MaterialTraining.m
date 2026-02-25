classdef MaterialTraining < handle
   
    properties (Access = public)
        designVariable
    end

    properties (Access = private)
        fHandle
        mesh
        mD
        inclusionType
        nSubdomains
        dim
        materialInterpolator
    end

    methods(Access=public)

        function obj=MaterialTraining(cParams)
            obj.init(cParams);
        end

        function material = create(obj)
            obj.mD = obj.createMeshDomain();

            switch obj.inclusionType
                case {'Hole','HoleRaul'}
                    [young,poisson] = obj.computeElasticProperties();
                    s.type          = 'ISOTROPIC';
                    s.ptype         = 'ELASTIC';
                    s.ndim          = obj.mesh.ndim;
                    s.young         = young;
                    s.poisson       = poisson;
                    material        = Material.create(s);
                case 'Material'
                    obj.createDesignVariable()
                    obj.createMaterialInterpolator(); 
                    s.type                 = 'DensityBased';
                    s.density              = obj.designVariable;
                    s.materialInterpolator = obj.materialInterpolator;
                    s.dim                  = obj.dim;
                    s.mesh                 = obj.mesh;
                    material = Material.create(s);
                    material.setDesignVariable({obj.designVariable.fun})
                    material = material.obtainTensor();
            end
        end

        function V=computeVolume(obj)
           V  = Integrator.compute(obj.designVariable.fun,obj.mesh,2);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access=private)
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.inclusionType  = cParams.inclusionType;
            obj.nSubdomains    = cParams.nSubdomains;
            obj.fHandle        = cParams.geomFun.getHandle;
            
            if obj.mesh.ndim == 2
                obj.dim='2D';
            elseif obj.mesh.ndim == 3
                obj.dim='3D';
            end

        end

        function [young,poisson] = computeElasticProperties(obj)
            E  = 1;
            nu  = 1/3;       
            young   = ConstantFunction.create(E,obj.mD);
            poisson = ConstantFunction.create(nu,obj.mD);
        end
      

        function meshDom = createMeshDomain(obj)
            if sum(obj.nSubdomains > 1)>= 1
                s.nsubdomains   = obj.nSubdomains; %nx ny
                s.meshReference = obj.mesh;
                s.tolSameNode   = 1e-10;
                m = MeshCreatorFromRVE2D(s);
                [meshDom,~,~,~,~,~,~] = m.create();
            else
                meshDom = obj.mesh;
            end
        end


        function createDesignVariable(obj)
            ls= obj.computeLevelSet();
            sUm.backgroundMesh = obj.mD;
            sUm.boundaryMesh   = obj.mD.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(-ls);
            funLS        = CharacteristicFunction.create(uMesh);
            s.filterType = 'LUMP';
            s.mesh       = obj.mD;
            s.trial      = LagrangianFunction.create(obj.mD,1,'P1');
            f = Filter.create(s);
            s.fun      = f.compute(funLS,2);
            s.type     = 'Density';
            s.plotting = false;
            dens = DesignVariable.create(s);
            obj.designVariable = dens;

        end

        function ls=computeLevelSet(obj)
            s.type   = 'Periodic';
            s.xmin   = min(obj.mesh.coord(:,1));
            s.xmax   = max(obj.mesh.coord(:,1));
            s.ymin   = min(obj.mesh.coord(:,2));
            s.ymax   = max(obj.mesh.coord(:,2));
            s.fBase  = obj.fHandle;
            g        = GeometricalFunction(s);
            phiFun   = g.computeLevelSetFunction(obj.mD);
            ls       = phiFun.fValues;
        end


        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            E1 = 1;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = obj.dim;
            s.matA           = matA;
            s.matB           = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator= m;
        end

    end

end

  




