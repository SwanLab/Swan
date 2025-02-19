classdef ConnectivityComputer < handle

    properties (Access = public)

    end

    properties (Access = private)
        mesh
        levelSet
        characteristicFunction
        designVariable

        materialInterpolation
        filter
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = ConnectivityComputer()
            obj.init();
            obj.createMesh();
            obj.createLevelSet();
            obj.createFilter();
            obj.createCharacteristicFunction();
            obj.createDesignVariable();

            % obj.levelSet.getUnfittedMesh().plot()
            % obj.density.plot()
            % shading flat
            % colormap('gray');
            % colormap(flipud(gray));
            % colorbar
            % figure

            obj.computeComplianceFunctional();
            obj.computeEigenValueFunctional()
        end

    end

    methods (Access = private)

        function init(obj)

        end

        function createMesh(obj)
            x1 = linspace(0,2,100);
            x2 = linspace(0,1,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
        end

        function createLevelSet(obj)
            s.type        = 'Circle';
            s.radius      = 0.3;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
            g             = GeometricalFunction(s);
            phi           = g.computeLevelSetFunction(obj.mesh);
            obj.levelSet = phi;
        end

        function createCharacteristicFunction(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(s);
            uMesh.compute(obj.levelSet.fValues);

            sC.uMesh = uMesh;
            obj.characteristicFunction  = CharacteristicFunction.create(sC);
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = P1Function.create(obj.mesh,1);
            f = Filter.create(s);
            obj.filter = f;
        end

        function createDesignVariable(obj)
            sD.fun  = obj.filter.compute(obj.characteristicFunction,'QUADRATIC');
            sD.mesh = obj.mesh;
            sD.type = 'Density';
            dens    = DesignVariable.create(sD);
            obj.designVariable = dens;
        end

        function computeComplianceFunctional(obj)
            obj.materialInterpolation = obj.computeMaterialInterpolation();
            s.mesh = obj.mesh;
            s.type    = 'ELASTIC';
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.bc = obj.createBoundaryConditionsForElasticity();
            s.interpolationType = 'LINEAR';
            fem = FEM.create(s);
            fem.solve();




            sH.type = 'ByInterpolation';
            sH.material = obj.createMaterial();
            sH.interpolationSettings.interpolation = obj.materialInterpolation;


            homogVarComp = HomogenizedVarComputer.create(sH);


            s.type = 'compliance';
            s.designVariable = obj.designVariable;
            s.femSettings.physicalProblem      = fem;
            s.femSettings.designVariableFilter = obj.filter;
            s.femSettings.gradientFilter       = obj.filter;
            s.homogVarComputer = homogVarComp;
            s.targetParameters = [];
            sh = ShapeFunctional.create(s);


            sh.computeFunctionAndGradient();
        end


        function matInt = computeMaterialInterpolation(obj)
            c.typeOfMaterial = 'ISOTROPIC';
            c.interpolation = 'SIMPALL';
            c.nElem = obj.mesh.nelem;
            c.dim = '2D';
            c.constitutiveProperties.rho_plus = 1;
            c.constitutiveProperties.rho_minus = 0;
            c.constitutiveProperties.E_plus = 1;
            c.constitutiveProperties.E_minus = 1e-3;
            c.constitutiveProperties.nu_plus = 1/3;
            c.constitutiveProperties.nu_minus = 1/3;

            matInt = MaterialInterpolation.create(c);
        end

        function mat = createMaterial(obj)
            d = obj.designVariable.fun.project('P0');
            %  d.fValues = round(d.fValues);

            dens  = d.fValues;
            mat  = obj.materialInterpolation.computeMatProp(dens);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = mat.kappa;
            s.mu    = mat.mu;
            mat = Material.create(s);
            mat.compute(s);
        end

        function bc = createBoundaryConditionsForElasticity(obj)
            bM = obj.mesh.createBoundaryMesh();

            dirichletBc.boundaryId=1;
            dirichletBc.dof=[1,2];
            dirichletBc.value=[0,0];
            newmanBc.boundaryId=2;
            newmanBc.dof=[2];
            newmanBc.value=[-1];

            [dirichlet,pointload] = obj.createBc(bM,dirichletBc,newmanBc);
            bc.dirichlet=dirichlet;
            bc.pointload=pointload;
        end

        function [dirichlet,pointload] = createBc(obj,boundaryMesh,dirchletBc,newmanBc)
            dirichlet = obj.createBondaryCondition(boundaryMesh,dirchletBc);
            pointload = obj.createBondaryCondition(boundaryMesh,newmanBc);
        end

        function cond = createBondaryCondition(obj,bM,condition)
            nbound = length(condition.boundaryId);
            cond = zeros(1,3);
            for ibound=1:nbound
                ncond  = length(condition.dof(nbound,:));
                nodeId= unique(bM{condition.boundaryId(ibound)}.globalConnec);
                nbd   = length(nodeId);
                for icond=1:ncond
                    bdcond= [nodeId, repmat(condition.dof(icond),[nbd,1]), repmat(condition.value(icond),[nbd,1])];
                    cond=[cond;bdcond];
                end
            end
            cond = cond(2:end,:);
        end

        function computeEigenValueFunctional(obj)
            eigen = obj.computeEigenValueProblem();
            s.eigenModes = eigen;
            s.designVariable = obj.designVariable;
            mE = MinimumEigenValueFunctional(s);
            mE.computeFunctionAndGradient()
        end

        function eigen = computeEigenValueProblem(obj)
            s.mesh = obj.mesh;
            eigen  = StiffnessEigenModesComputer(s);
        end



    end

end