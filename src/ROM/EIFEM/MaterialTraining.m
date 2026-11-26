classdef MaterialTraining < handle
   
    properties (Access = public)
        designVariable
        levelSet
    end

    properties (Access = private)
        fHandle
        mesh
        dim
        bounds
        E1
        nu1
        E0
        nu0
    end

    methods(Access=public)

        function obj=MaterialTraining(cParams)
            obj.init(cParams);
        end

        function material = create(obj)
            obj.createDesignVariable()
            mat = obj.createMaterialInterpolator();
            s.type                 = 'DensityBased';
            s.density              = obj.designVariable;
            material = mat(s.density.fun);
        end

        function V=computeVolume(obj)
            V  = Integrator.compute(obj.designVariable.fun,obj.mesh,2);
        end
    end

    methods(Access=private)
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.fHandle        = cParams.geomFun.getHandle;
            obj.bounds         = cParams.bounds;
            obj.E1             = cParams.E;
            obj.nu1            = cParams.nu;
            obj.E0             = 1e-3;
            obj.nu0            = 1/3;
        end

        function createDesignVariable(obj)
            ls= obj.computeLevelSet();

            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);

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
                    s.zmin   = 0;
                    s.zmax   = obj.bounds(3);
            end
            s.xmin   = 0;
            s.xmax   = obj.bounds(1);
            s.ymin   = 0;
            s.ymax   = obj.bounds(2);
            s.fBase  = obj.fHandle;
            g        = GeometricalFunction(s);
            phiFun   = g.computeLevelSetFunction(obj.mesh);
            obj.levelSet = phiFun;
            ls       = phiFun.fValues;
        end


        function C=createMaterialInterpolator(obj)
            N      = obj.mesh.ndim;
            muA    = obj.computeMu(obj.E0,obj.nu0);
            kappaA = obj.computeKappa(obj.E0,obj.nu0,N);
            muB    = obj.computeMu(obj.E1,obj.nu1);
            kappaB = obj.computeKappa(obj.E1,obj.nu1,N);

            mu     = @(rho) SimpAllInterpolator.computeMu(muA,muB,kappaA,kappaB,rho,N); mu = @(rho) Expand(mu(rho),4);
            kappa  = @(rho) SimpAllInterpolator.computeKappa(muA,muB,kappaA,kappaB,rho,N); kappa = @(rho) Expand(kappa(rho),4);
            lambda = @(rho) kappa(rho) - (2/N)*mu(rho);
            I      = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI    = ConstantFunction.create(kronEye(N),obj.mesh);
            C      = @(rho) 2*mu(rho).*I + lambda(rho).*IxI;
        end

        function mu = computeMu(~,E,nu)
            mu = E./(2*(1+nu));
        end

        function kappa = computeKappa(~,E,nu,N)
            kappa = E./(N*(1-(N-1)*nu));
        end

    end

end
