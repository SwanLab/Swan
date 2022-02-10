classdef NewElasticProblem < handle %NewFEM

    properties (GetAccess = public, SetAccess = private)
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
        % Kgen %"LHS"
    end

    properties (Access = private)
        fileName
        femData
        mesh
        problemData
        dof
        stiffnessMatrix
        forces
        solver
    end

    methods (Access = public)

        function obj = NewElasticProblem(cParams)
            obj.init(cParams);
            obj.readProblemData();
            obj.createQuadrature();
            obj.computeDimensions();
            obj.createMaterial();
            obj.computeMaterialProperties();
            obj.createInterpolation();
            obj.createBCApplier();
            obj.createSolver();
        end

        
        function computeVariables(obj)
            obj.computeStiffnessMatrix()
            obj.computeForces();
            obj.computeDisplacements();
        end
    
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.fileName = cParams.fileName;
            obj.nFields = 1;
        end

        function readProblemData(obj)
            obj.createFemData();
            obj.createProblemData();
            obj.mesh = obj.femData.mesh;
        end

        function createFemData(obj)
            fName       = obj.fileName;
            femReader   = FemInputReader_GiD();
            obj.femData = femReader.read(fName);
        end

        function createProblemData(obj)
            s = obj.femData;
            pd.fileName     = obj.fileName;
            pd.scale        = s.scale;
            pd.pdim         = s.pdim;
            pd.ptype        = s.ptype;
            pd.nelem        = s.mesh.nelem;
            pd.bc.dirichlet = s.dirichlet;
            pd.bc.pointload = s.pointload;
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
            %int.computeShapeDeriv(obj.quadrature.posgp);
            obj.interp{1} = int;
        end

        function createBCApplier(obj)
            obj.dof = DOF_Elastic(obj.fileName,obj.mesh,obj.problemData.pdim,obj.nFields,obj.interp);
            cParams.nfields = obj.nFields;
            cParams.dof     = obj.dof;
            cParams.scale   = obj.problemData.scale;
            cParams.type    = 'Dirichlet'; % defined in Element
            obj.bcApplier = BoundaryConditionsApplier.create(cParams);
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
            obj.stiffnessMatrix = Kred;
        end

        function computeForces(obj)
            f    = obj.computeExternalForces();
            fRed = obj.reduceForcesMatrix(f);
            obj.forces = fRed;
%             R = obj.compute_imposed_displacement_force(obj.K);
%             obj.fext = Fext + R;
%             obj.rhs = obj.integrator.integrate(fNodal);
        end

        function f = computeExternalForces(obj)
            s.dim  = obj.dim;
            s.dof  = obj.dof;
            s.mesh = obj.mesh;
            fcomp = ForcesComputer(s);
            f = fcomp.compute();
        end

        function Fred = reduceForcesMatrix(obj, forces)
            Fred = obj.bcApplier.fullToReducedVector(forces);
        end
       
        function u = computeDisplacements(obj)
            Kred = obj.stiffnessMatrix;
            Fred = obj.forces;
            u = obj.solver.solve(Kred,Fred);
            u = obj.bcApplier.reducedToFullVector(u);
            obj.variables.d_u = u;
        end

    end

end