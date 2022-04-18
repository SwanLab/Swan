classdef DiffReactProblem < handle
    
    properties (GetAccess = public, SetAccess = protected)
        variables
    end
    
    properties (Access = protected)
        dim
        M,K
        mesh
        solver
        epsilon
        problemData
        boundaryConditions
    end

    methods (Static, Access = public)

        function obj = create(s)
            isRobin = isfield(s, 'isRobinTermAdded') ...
                      && s.isRobinTermAdded == 1;
            switch isRobin
                case true
                    obj = DiffReactProblem_Robin(s);
                case false
                    obj = DiffReactProblem_Neumann(s);
            end
        end

    end

    methods (Access = public)
        
        function obj = DiffReactProblem(cParams)
            obj.init(cParams);
            obj.computeProblemDimensions();
            obj.createBoundaryConditions();
            obj.createSolver();
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix();
        end

        function obj = setEpsilon(obj,epsilon)
            obj.epsilon = epsilon;
        end

        function M = getM(obj)
            M = obj.M;
        end

        function K = getK(obj)
            K = obj.K;
        end

        function dvol = computeDvolume(obj)
            int = obj.mesh.interpolation;
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            dvol = g.dvolu;
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh = cParams.mesh;
            obj.problemData.pdim = '1D';
            obj.problemData.scale = cParams.scale;
            if isfield(cParams,'fileName') % Robin
                obj.problemData.fileName = cParams.fileName;
            end
        end
        
        function computeProblemDimensions(obj)
            m = obj.mesh;
            obj.dim = obj.computeDimensions(m);
        end

        function d = computeDimensions(obj, mesh)
            s.ngaus = [];
            s.mesh  = mesh;
            s.pdim  = obj.problemData.pdim;
            dims    = DimensionVariables(s);
            dims.compute();
            d = dims;
        end

        function createBoundaryConditions(obj)
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.scale        = obj.problemData.scale;
            s.bc.dirichlet = [];
            s.bc.pointload = [];
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
            obj.K = LHS.compute();
        end
        
        function computeMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end
    
    end
    
end