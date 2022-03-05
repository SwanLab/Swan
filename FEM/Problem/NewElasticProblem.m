classdef NewElasticProblem < handle %NewFEM

%     properties (GetAccess = public, SetAccess = private)
    properties (Access = public)
        variables
    end

    properties (Access = private)
        material
        nFields
        interp
        bcApplier

        % Poda 1
        quadrature
        dim
        boundaryConditions
        displacement
    end

    properties (Access = private)
        mesh
        problemData
        stiffnessMatrix
        stiffnessMatrixRed
        forces
        solver
    end

    methods (Access = public)

        function obj = NewElasticProblem(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.computeDimensions();
            obj.createMaterial();
            obj.computeMaterialProperties();
            obj.createInterpolation();
            obj.createBoundaryConditions();
            obj.createBCApplier();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.computeDisplacements();
            obj.computeStrain();
            obj.computeStress();
            obj.computePrincipalDirection();
        end

        function plot(obj)
            s.dim            = obj.dim;
            s.mesh           = obj.mesh;
            s.displacement = obj.variables.d_u;
            plotter = FEMPlotter(s);
            plotter.plot();
        end

        function dim = getDimensions(obj)
            dim = obj.dim;
        end

        function setC(obj, C)
            obj.material.C = C;
        end

        function dvolu = getDvolume(obj)
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.nFields = 1;
            obj.mesh        = cParams.mesh;
            pd.scale        = cParams.scale;
            pd.pdim         = cParams.dim;
            pd.ptype        = cParams.type;
            pd.bc.dirichlet = cParams.dirichlet;
            pd.bc.pointload = cParams.pointload;
            obj.problemData = pd;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function computeDimensions(obj)
            s.ngaus = obj.quadrature.ngaus;
            s.mesh  = obj.mesh;
            s.pdim  = obj.problemData.pdim;
            d       = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end

        function createMaterial(obj)
            s.ptype = obj.problemData.ptype;
            s.pdim  = obj.problemData.pdim;
            s.nelem = obj.mesh.nelem;
            %s.geometry = obj.geometry; % Hyperelastic
            s.mesh  = obj.mesh;
            obj.material = Material.create(s);
        end

        function computeMaterialProperties(obj)
            I = ones(obj.dim.nelem,obj.dim.ngaus);
            s.kappa = .9107*I;
            s.mu    = .3446*I;
            obj.material.compute(s);
        end

        function createInterpolation(obj)
            int = obj.mesh.interpolation;
            obj.interp{1} = int;
        end

        function createBoundaryConditions(obj)
            s.dim          = obj.dim;
            s.bc           = obj.problemData.bc;
            s.globalConnec = obj.mesh.connec;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createBCApplier(obj)
            s.BC      = obj.boundaryConditions;
            s.dim     = obj.dim;
            s.nfields = obj.nFields;
            s.scale   = obj.problemData.scale;
            s.type    = 'Dirichlet'; % defined in Element
            obj.bcApplier = BoundaryConditionsApplier.create(s);
        end

        function createSolver(obj)
            obj.solver = Solver.create;
        end

        function computeStiffnessMatrix(obj)
            s.type = 'ElasticStiffnessMatrixOld';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            s.material     = obj.material;
            LHS = LHSintegrator.create(s);
            K   = LHS.compute();
            Kred = obj.bcApplier.fullToReducedMatrix(K);
            obj.stiffnessMatrix    = K;
            obj.stiffnessMatrixRed = Kred;
        end

        function computeForces(obj)
            f    = obj.computeExternalForces();
            fRed = obj.reduceForcesMatrix(f);
            obj.forces = fRed;
            %             R = obj.compute_imposed_displacement_force(obj.K);
            %             obj.fext = Fext + R;
            %             obj.rhs = obj.integrator.integrate(fNodal);
        end

        function F = computeExternalForces(obj)
            s.dim  = obj.dim;
            s.BC   = obj.boundaryConditions;
            s.mesh = obj.mesh;
            fcomp = ForcesComputer(s);
            F = fcomp.compute();
            R = fcomp.computeReactions(obj.stiffnessMatrix);
            obj.variables.fext = F + R;
        end

        function Fred = reduceForcesMatrix(obj, forces)
            Fred = obj.bcApplier.fullToReducedVector(forces);
        end

        function u = computeDisplacements(obj)
            Kred = obj.stiffnessMatrixRed;
            Fred = obj.forces;
            u = obj.solver.solve(Kred,Fred);
            u = obj.bcApplier.reducedToFullVector(u);
            obj.variables.d_u = u;
        end

        function computeStrain(obj)
            s.dim                = obj.dim;
            s.mesh               = obj.mesh;
            s.quadrature         = obj.quadrature;
            s.displacement       = obj.variables.d_u;
            s.interpolation      = obj.interp{1};
            s.boundaryConditions = obj.boundaryConditions;
            scomp  = StrainComputer(s);
            strain = scomp.compute();
            obj.variables.strain = strain;
        end

        function computeStress(obj)
            s.C      = obj.material.C;
            s.dim    = obj.dim;
            s.strain = obj.variables.strain;
            scomp  = StressComputer(s);
            stress = scomp.compute();
            obj.variables.stress = stress;
        end

        function computePrincipalDirection(obj)
            stress = obj.variables.stress;
            s.type             = obj.problemData.pdim;
            s.eigenValueComputer.type = 'PRECOMPUTED';
            pcomp = PrincipalDirectionComputer.create(s);
            pcomp.compute(stress);
            obj.variables.principalDirections = pcomp.direction;
            obj.variables.principalStress     = pcomp.principalStress;
        end

    end

end