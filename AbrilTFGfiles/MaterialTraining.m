classdef MaterialTraining < handle
   
    properties (Access = public)
        designVariable
        levelSet
        unfittedMesh
    end

    properties (Access = private)
        fHandle
        mesh
        dim
        materialInterpolator
    end

    methods(Access=public)

        function obj=MaterialTraining(cParams)
            obj.init(cParams);
        end

        function material = create(obj)
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

        function V=computeVolume(obj)
           V  = Integrator.compute(obj.designVariable.fun,obj.mesh,2);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access=private)
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.fHandle        = cParams.geomFun.getHandle;
            
            if obj.mesh.ndim == 2
                obj.dim='2D';
            elseif obj.mesh.ndim == 3
                obj.dim='3D';
            end

        end

        function createDesignVariable(obj)
            ls= obj.computeLevelSet();
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(-ls);
            obj.unfittedMesh = uMesh;
            Mprint=UnfittedMesh(sUm);
            Mprint.compute(ls);
            funLS        = CharacteristicFunction.create(uMesh);
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            s.fun      = f.compute(funLS,2);
            s.type     = 'Density';
            s.plotting = false;
            dens = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function ls=computeLevelSet(obj)
            switch obj.mesh.ndim
                case 2
                    s.type   = 'Periodic';
                case 3
                    s.type   = 'Periodic3D';
                    s.zmin   = -1;
                    s.zmax   = 1;
            end
            s.xmin   = -1;
            s.xmax   = 1;
            s.ymin   = -1;
            s.ymax   = 1;
            s.fBase  = obj.fHandle;
            g        = GeometricalFunction(s);
            phiFun   = g.computeLevelSetFunction(obj.mesh);
            obj.levelSet = phiFun;
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
