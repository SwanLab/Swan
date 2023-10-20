classdef PhaseFieldComputer < handle

    properties (Access = public)

    end

    properties (Access = private)
        mesh
        boundaryConditions
        materialInterpolation
        material
        phaseField
        fem
    end


    methods (Access = public)

        function obj = PhaseFieldComputer()
            obj.init()
            obj.computeMesh();
            obj.createBoundaryConditions();
            obj.createPhaseField();
            obj.createMaterialInterpolation();
            obj.computeFEM();
            obj.createFGaussEnergyFunction();
        end

    end




    methods (Access = private)

        function init(obj)

        end

        function computeMesh(obj)
            %file = 'test2d_triangle';
            %a.fileName = file;
            % s = FemDataContainer(a);

            % Generate coordinates
            x1 = linspace(0,2,20);
            x2 = linspace(1,2,20);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord = V(:,1:2);
            s.connec = F;
            m = Mesh(s);
            obj.mesh = m;
        end

        function createBoundaryConditions(obj)
            dirichletNodes = abs(obj.mesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(obj.mesh.coord(:,1));
            isInRight = abs(obj.mesh.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(obj.mesh.coord(:,2)-1.5) < 0.1;
            forceNodes = isInRight & isInMiddleEdge;
            nodes = 1:obj.mesh.nnodes;

            ndim = 2;
            bc.dirichlet = zeros(ndim*length(nodes(dirichletNodes)),3);
            for i=1:ndim
                bc.dirichlet(i:2:end,1) = nodes(dirichletNodes);
                bc.dirichlet(i:2:end,2) = i;
            end
            bc.dirichlet(:,3) = 0;

            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = -1;
            obj.boundaryConditions = bc;
        end

        function createPhaseField(obj)
            % sLS.type       = 'rectangleInclusion';
            % sLS.mesh       = obj.mesh;
            % sLS.ndim       = 2;
            % sLS.widthH = 1;
            % sLS.widthV = 2;
            % sLS.coord      = obj.mesh.coord;
            % ls = LevelSetCreator.create(sLS);
            % dmg = ls.getValue();
            % obj.phi = dmg;
            % scatter3(obj.mesh.coord(:,1),obj.mesh.coord(:,2),obj.phi);

            xmax = max(obj.mesh.coord(:,1));
            sAF.fHandle = @(x) x(1,:,:)/xmax;
            sAF.ndimf   = 1;
            sAF.mesh    = obj.mesh;
            xFun = AnalyticalFunction(sAF);

            phi = xFun.project('P1');
            phi.plot();
            obj.phaseField = phi;
        end

        function createMaterialInterpolation(obj)
            c.typeOfMaterial = 'ISOTROPIC';
            c.interpolation = 'PhaseFieldI';
            c.nElem = obj.mesh.nelem;
            c.dim = '2D';
            c.constitutiveProperties.rho_plus = 1;
            c.constitutiveProperties.rho_minus = 0;
            c.constitutiveProperties.E_plus = 1;
            c.constitutiveProperties.E_minus = 1e-3;
            c.constitutiveProperties.nu_plus = 1/3;
            c.constitutiveProperties.nu_minus = 1/3;

            matInt = MaterialInterpolation.create(c);
            obj.materialInterpolation = matInt;
        end


        function mat = createMaterial(obj)
            phi0 = obj.phaseField.project('P0');
            mat  = obj.materialInterpolation.computeMatProp(squeeze(phi0.fValues));
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = mat.kappa;
            s.mu    = mat.mu;
            mat = Material.create(s);
            mat.compute(s);
            obj.material = mat;
        end


        function computeFEM(obj)
            s.mesh = obj.mesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.bc = obj.boundaryConditions;
            obj.fem = FEM.create(s);
            obj.fem.solve();

            obj.fem.uFun.plot()
            obj.fem.stressFun.plot()
            obj.fem.strainFun.plot()

            obj.fem.uFun.fValues(:,end+1) = 0;
            obj.fem.uFun.ndimf = 3;
            obj.fem.print('Example','Paraview')
        end

        function createFGaussEnergyFunction(obj)
            e = obj.fem.strainFun;
            s.fValues = sum(e.fValues.*e.fValues);
            s.quadrature = e.quadrature;
            s.mesh = obj.mesh;
            energy = FGaussDiscontinuousFunction(s);

            energy.plot();
        end

        function createEnergyMassMatrix()
            s.function =  obj.createFGaussEnergyFunction()

            LHS
            M = LHS.create(); 
        end


    end

end
