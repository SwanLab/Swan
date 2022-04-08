classdef ThermalProblem < handle
    
    properties (Access = public)
        variables
    end

    properties (Access = private)
        nFields
        boundaryConditions
        displacement
        problemData
        stiffnessMatrix
        stiffnessMatrixRed
        forces
        solver
        geometry

        dofsInElem
    end

    properties (Access = protected)
        quadrature
        dim
        material

        vstrain

        mesh, interp % For Homogenization
    end

    methods (Access = public)

        function obj = ThermalProblem(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.createQuadrature();
            obj.createGeometry();
            obj.computeDimensions();
            obj.createBoundaryConditions();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.computeDisplacements();
        end

        function plot(obj)
            s.dim            = obj.dim;
            s.mesh           = obj.mesh;
            s.displacement = obj.variables.d_u;
            plotter = FEMPlotter(s);
            plotter.plot();
        end

        function dvolu = getDvolume(obj)
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
        end

        function print(obj,fileName)
            dI = obj.createPostProcessDataBase(fileName);
            dI.pdim = '2D';
            dI.ptype = 'THERMAL';
            dI.name = 'temp';
            postprocess = Postprocess('ScalarNodal',dI);
            q = obj.quadrature;
            d.fields = obj.variables.d_u;
            d.quad = q;
            postprocess.print(1,d);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.nFields = 1;
            obj.mesh        = cParams.mesh;
            pd.scale        = cParams.scale;
            pd.pdim         = '1D';
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

        function createInterpolation(obj)
            int = obj.mesh.interpolation;
            obj.interp{1} = int;
        end
       
        function createGeometry(obj)
            q   = obj.quadrature;
            int = obj.interp{1};
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

        function createBoundaryConditions(obj)
            s.dim        = obj.dim;
            s.mesh       = obj.mesh;
            s.scale      = obj.problemData.scale;
            s.bc         = obj.problemData.bc;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createSolver(obj)
            obj.solver = Solver.create();
        end

        function computeStiffnessMatrix(obj)
            s.type = 'StiffnessMatrix';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            K   = LHS.compute();
            Kred = obj.boundaryConditions.fullToReducedMatrix(K);
            obj.stiffnessMatrix    = K;
            obj.stiffnessMatrixRed = Kred;
        end

        function computeForces(obj)
            f    = obj.computeExternalForces();
            fRed = obj.reduceForcesMatrix(f);
            obj.forces = fRed;
        end

        function F = computeExternalForces(obj)
            s.dim         = obj.dim;
            s.BC          = obj.boundaryConditions; %dofsInElem, Neumann
            s.mesh        = obj.mesh;
            s.material    = obj.material;
            s.geometry    = obj.geometry;
            s.dvolume     = obj.getDvolume();
            s.globalConnec = obj.mesh.connec;
            if isprop(obj, 'vstrain')
                s.vstrain = obj.vstrain;
            end
            fcomp = ForcesComputer(s);
            F = fcomp.compute();
            R = fcomp.computeReactions(obj.stiffnessMatrix);
            obj.variables.fext = F + R;
        end

        function Fred = reduceForcesMatrix(obj, forces)
            Fred = obj.boundaryConditions.fullToReducedVector(forces);
        end

        function u = computeDisplacements(obj)
            Kred = obj.stiffnessMatrixRed;
            Fred = obj.forces;
            u = obj.solver.solve(Kred,Fred);
            u = obj.boundaryConditions.reducedToFullVector(u);
            obj.variables.d_u = u;
        end

        function d = createPostProcessDataBase(obj,fileName)
            dI.mesh    = obj.mesh;
            dI.outName = fileName;
            dI.pdim = '2D';
            dI.ptype = 'THERMAL';
            ps = PostProcessDataBaseCreator(dI);
            d = ps.getValue();
        end
        
    end

end